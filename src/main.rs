use std::{
    iter::Sum,
    ops::{Add, Mul, Sub},
};

use fraction::{GenericFraction, ToPrimitive};
use indicatif::{ProgressIterator, ProgressStyle};
use num::BigInt;
use rand::Rng;
use rayon::prelude::*;

use plotters::{prelude::*, style::full_palette::PINK};

type Error = Box<dyn std::error::Error>;
type Fraction = GenericFraction<BigInt>;
type SmallFraction = GenericFraction<u64>;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
struct Position {
    x: isize,
    y: isize,
    z: isize,
}

#[derive(Debug, PartialEq, Clone)]
struct FPosition {
    x: SmallFraction,
    y: SmallFraction,
    z: SmallFraction,
}

impl From<Position> for FPosition {
    fn from(value: Position) -> Self {
        FPosition {
            x: value.x.into(),
            y: value.y.into(),
            z: value.z.into(),
        }
    }
}

impl From<(isize, isize, isize)> for Position {
    fn from((x, y, z): (isize, isize, isize)) -> Self {
        Position { x, y, z }
    }
}

impl Position {
    fn sqr_len(&self) -> isize {
        let squared_x = self.x * self.x;
        let squared_y = self.y * self.y;
        let squared_z = self.z * self.z;
        squared_x + squared_y + squared_z
    }
}

impl Mul<FPosition> for SmallFraction {
    type Output = FPosition;

    fn mul(self, rhs: FPosition) -> Self::Output {
        FPosition {
            x: self * rhs.x,
            y: self * rhs.y,
            z: self * rhs.z,
        }
    }
}

impl Sum for Position {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Position { x: 0, y: 0, z: 0 }, |acc, pos| acc + pos)
    }
}

impl Add for Position {
    type Output = Position;

    fn add(self, rhs: Self) -> Self::Output {
        Position {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            z: self.z + rhs.z,
        }
    }
}

#[derive(Debug)]
struct Lattice {
    sites_x: usize,
    sites_y: usize,
    sites_z: usize,
    start_x: usize,
    start_y: usize,
    start_z: usize,
}

impl Lattice {
    fn new(x_size: usize, y_size: usize, z_size: usize) -> Self {
        Lattice {
            sites_x: x_size - 1,
            sites_y: y_size - 1,
            sites_z: z_size - 1,
            start_x: x_size / 2,
            start_y: y_size / 2,
            start_z: z_size / 2,
        }
    }

    fn starting_point(&self) -> Position {
        (
            self.start_x as isize,
            self.start_y as isize,
            self.start_z as isize,
        )
            .into()
    }

    fn get_size(&self) -> (usize, usize, usize) {
        (self.sites_x, self.sites_y, self.sites_z)
    }

    fn adjacent(&self, &Position { x, y, z }: &Position) -> Vec<Position> {
        let mut avec = Vec::with_capacity(4);
        if x < self.sites_x as isize {
            // not on the right edge
            avec.push(Position { x: x + 1, y, z });
        };
        if x > 0 {
            // not on the left edge
            avec.push(Position { x: x - 1, y, z });
        };
        if y < self.sites_y as isize {
            // not on the bottom edge
            avec.push(Position { x, y: y + 1, z });
        };
        if y > 0 {
            // not on the top edge
            avec.push(Position { x, y: y - 1, z });
        };
        if z < self.sites_z as isize {
            // not on the top
            avec.push(Position { x, y, z: z + 1 });
        };
        if z > 0 {
            // not on the bottom
            avec.push(Position { x, y, z: z - 1 });
        };
        avec
    }
}

impl Sub for FPosition {
    type Output = FPosition;

    fn sub(self, rhs: Self) -> Self::Output {
        FPosition {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

impl Sub for Position {
    type Output = Position;

    fn sub(self, rhs: Self) -> Self::Output {
        Position {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            z: self.z - rhs.z,
        }
    }
}

#[derive(Debug)]
struct Polymer<'a> {
    positions: Vec<Position>,
    weights: Vec<Fraction>,
    head: Position,
    lattice: &'a Lattice,
    len: usize,
}

impl Polymer<'_> {
    fn weight(&self, length: usize) -> &Fraction {
        &self.weights[length]
    }
}

enum GrowthResult {
    Success,
    Cornered,
}

trait Observable {
    fn value(polymer: &Polymer, length: usize) -> Fraction;

    fn weighted(polymers: &[Polymer], length: usize, total_weight: Fraction) -> Fraction {
        let numerator: Fraction = polymers
            .par_iter()
            .filter(|poly| poly.len >= length)
            .map(|poly| &Self::value(poly, length) * poly.weight(length))
            .sum();
        numerator / total_weight
    }
}

struct Gyration;
impl Gyration {
    fn pdif_center(
        &Position {
            x: px,
            y: py,
            z: pz,
        }: &Position,
        &FPosition {
            x: cx,
            y: cy,
            z: cz,
        }: &FPosition,
    ) -> Fraction {
        let diff_x = Into::<SmallFraction>::into(px) - cx;
        let diff_y = Into::<SmallFraction>::into(py) - cy;
        let diff_z = Into::<SmallFraction>::into(pz) - cz;
        let sdiff = diff_x * diff_x + diff_y * diff_y + diff_z * diff_z;

        sdiff.into_fraction()
    }
}

impl Observable for Gyration {
    fn value(polymer: &Polymer, length: usize) -> Fraction {
        let pref =
            Into::<SmallFraction>::into(1) / Into::<SmallFraction>::into(polymer.positions.len());
        let center: FPosition = pref
            * FPosition::from(
                polymer
                    .positions
                    .iter()
                    .take(length)
                    .map(|&pos| pos)
                    .sum::<Position>(),
            );
        pref.into_fraction::<BigInt>()
            * polymer
                .positions
                .iter()
                .take(length)
                .map(|particle| Self::pdif_center(particle, &center))
                .sum::<Fraction>()
    }
}

struct EndToEnd;
impl Observable for EndToEnd {
    fn value(polymer: &Polymer, length: usize) -> Fraction {
        (polymer.positions[length] - polymer.positions[0])
            .sqr_len()
            .into()
    }
}

impl<'a> Polymer<'a> {
    fn grow(&mut self) -> GrowthResult {
        let adj = self.lattice.adjacent(&self.head);
        let valid_adj: Vec<_> = adj
            .iter()
            .filter(|pos| !self.positions.contains(pos))
            .collect();
        // if all surrounding positions are in the body we can't move
        if valid_adj.len() == 0 {
            return GrowthResult::Cornered;
        }
        self.len += 1;
        self.weights
            .push(self.weights.last().expect("starting weight exists") * valid_adj.len());
        // pick a random next position
        let next = valid_adj[rand::rng().random_range(0..valid_adj.len())];
        // set the new head
        self.head = *next;
        // add the new head to the body
        self.positions.push(self.head);

        GrowthResult::Success
    }

    fn new(lattice: &'a Lattice) -> Polymer<'a> {
        let start = lattice.starting_point();
        Polymer {
            positions: vec![start],
            head: start,
            weights: vec![1.into()],
            lattice,
            len: 0,
        }
    }

    fn grow_perm_till_end(
        lattice: &'a Lattice,
        w_low: usize,
        w_high: usize,
        max_len: usize,
    ) -> Vec<Polymer<'a>> {
        let mut polies = vec![Polymer::new(&lattice)];

        polies.iter().flat_map(|poly| {
            if poly.len < w_low {
                if rand::random::<f32>() < 0.5 {
                    vec![]
                } else {
                    poly.weights = poly.weights * 2.into();
                    vec![poly]
                }
            } else if poly.len > w_high {
                poly.weights = poly.weights / 2.into();
                vec![poly, poly]
            } else {
                vec![poly]
            }
        });

        polies
    }

    fn grow_till_end(&mut self) {
        while let GrowthResult::Success = self.grow() {}
        if self.len < W_LOWER {
            if rand::random::<f32>() < 0.5 {
                self.state = PERMState::SmallDitch
            } else {
                self.state = PERMState::SmallDoubleKeep
            }
        } else if self.len > W_UPPER {
            self.state = PERMState::Large
        } else {
            self.state = PERMState::Medium
        }
    }
}

enum Waa<'a> {
    Lower,
    Higher(Polymer<'a>),
    Between,
}

const WALK_AMOUNT: usize = 10;
const LATTICE_SIZE: usize = 350;

const W_LOWER: usize = 10;
const W_UPPER: usize = 40;

fn plot_observables(
    e2es: Vec<Fraction>,
    gyros: Vec<Fraction>,
    amounts: Vec<usize>,
    max_len: usize,
) -> Result<(), Error> {
    let root = BitMapBackend::new("observe.png", (1024, 768)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption("Observables", ("sans-serif", 50).into_font())
        .margin(5)
        .set_all_label_area_size(50)
        .build_cartesian_2d(
            0..max_len,
            0.0..gyros
                .iter()
                .map(|gyro| ToPrimitive::to_f64(gyro).unwrap_or(f64::NAN))
                .chain(
                    e2es.iter()
                        .map(|e2e| ToPrimitive::to_f64(e2e).unwrap_or(f64::NAN)),
                )
                // this max isn't fully correct I think
                .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                .expect("there should be a gyro or e2e"),
        )?;

    chart.configure_mesh().draw()?;

    chart.draw_series(LineSeries::new(
        (0..max_len).map(|x| (x, ToPrimitive::to_f64(&e2es[x]).expect("waa"))),
        RED.stroke_width(3),
    ))?;
    chart.draw_series(LineSeries::new(
        (0..max_len).map(|x| (x, ToPrimitive::to_f64(&gyros[x]).expect("woo"))),
        BLUE.stroke_width(3),
    ))?;
    chart.draw_series(LineSeries::new(
        (0..max_len).map(|x| (x, amounts[x] as f64)),
        PINK.stroke_width(3),
    ))?;
    root.present()?;

    Ok(())
}

fn plot_walks(size: (usize, usize, usize), polies: &[Polymer]) -> Result<(), Error> {
    // create a plot with all of them drawn
    let root = BitMapBackend::new("walk.png", (640, 640)).into_drawing_area();
    root.fill(&WHITE)?;
    let mut chart = ChartBuilder::on(&root)
        .caption(
            "Random non self-intersecting walk",
            ("sans-serif", 50).into_font(),
        )
        .margin(5)
        .set_all_label_area_size(50)
        .build_cartesian_3d(0..size.0, 0..size.1, 0..size.2)?;

    chart.configure_axes().draw()?;
    // draw each polymer
    for walk in (0..WALK_AMOUNT).progress_with_style(ProgressStyle::default_bar().template("Drawing plots [{elapsed_precise}/{duration_precise} (Remaining: {eta})] {bar:40.cyan/blue} {pos:>7}/{len:7} {per_sec} {msg}")?) {
        let part = walk as f32 / WALK_AMOUNT as f32;
        let color = ViridisRGBA::get_color(part);
        chart.draw_series(LineSeries::new(
            polies[walk]
                .positions
                .iter()
                .map(|&Position { x, y, z }| (x as usize, y as usize, z as usize)),
            color.stroke_width(3),
        ))?;
    }
    root.present()?;

    Ok(())
}

fn main() -> Result<(), Error> {
    // initialize the lattice
    let lattice = Lattice::new(LATTICE_SIZE, LATTICE_SIZE, LATTICE_SIZE);

    // create and walk all of the polymers
    let mut polies: Vec<Polymer> = (0..WALK_AMOUNT).map(|_| Polymer::new(&lattice)).collect();
    polies.iter_mut().progress_with_style(ProgressStyle::default_bar().template("Generating polymers [{elapsed_precise}/{duration_precise} (Remaining: {eta})] {bar:40.cyan/blue} {pos:>7}/{len:7} {per_sec} {msg}")?).par_bridge().for_each(|poly| {
        poly.grow_till_end();
    });

    plot_walks(lattice.get_size(), &polies)?;

    let max_poly_len = polies
        .par_iter()
        .max_by(|a, b| a.len.cmp(&b.len))
        .expect("there is a polymer")
        .len;
    let sum_weights: Vec<Fraction> = (0..max_poly_len)
        .progress_with_style(ProgressStyle::default_bar().template("Weights [{elapsed_precise}/{duration_precise}] ({eta} left) {barr:40.cyan/blue} {pos:>7}/{len:7} {per_sec} {msg}")?)
        .map(|length| {
            polies
                .par_iter()
                .filter(|poly| poly.len >= length)
                .map(|poly| poly.weight(length).clone())
                .sum()
        })
        .collect();

    let e2es: Vec<_> = (0..max_poly_len)
       .progress_with_style(ProgressStyle::default_bar().template("End to End [{elapsed_precise}/{duration_precise} (Remaining: {eta})] {bar:40.cyan/blue} {pos:>7}/{len:7} {per_sec} {msg}")?) 
        .map(|len| EndToEnd::weighted(&polies, len, sum_weights[len].clone()))
        .collect();
    let gyros: Vec<_> = (0..max_poly_len)
       .progress_with_style(ProgressStyle::default_bar().template("Gyration [{elapsed_precise}/{duration_precise} (Remaining: {eta})] {bar:40.cyan/blue} {pos:>7}/{len:7} {per_sec} {msg}")?) 
        .map(|len| Gyration::weighted(&polies, len, sum_weights[len].clone()))
        .collect();
    let amounts: Vec<_> = (0..max_poly_len)
        .map(|len| {
            polies
                .iter()
                .filter(|poly| poly.len >= len)
                .enumerate()
                .last()
                .unwrap()
                .0
        })
        .collect();

    for (idx, (gyro, e2e)) in gyros.iter().zip(e2es.iter()).enumerate() {
        println!("{} gyro: {:.2}, e2e: {:.2}", idx, gyro, e2e);
    }

    plot_observables(e2es, gyros, amounts, max_poly_len)?;

    Ok(())
}
