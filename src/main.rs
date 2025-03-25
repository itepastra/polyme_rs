use std::{
    iter::Sum,
    ops::{Add, Mul, Sub},
};

use fraction::{GenericFraction, ToPrimitive};
use indicatif::{ProgressBar, ProgressIterator};
use num::{BigInt, BigUint};
use rand::Rng;
use rayon::prelude::*;

use plotters::prelude::*;

type Error = Box<dyn std::error::Error>;
type Fraction = GenericFraction<BigInt>;
type SmallFraction = GenericFraction<u64>;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
struct Position {
    x: isize,
    y: isize,
}

#[derive(Debug, PartialEq, Clone)]
struct FPosition {
    x: SmallFraction,
    y: SmallFraction,
}

impl From<Position> for FPosition {
    fn from(value: Position) -> Self {
        FPosition {
            x: value.x.into(),
            y: value.y.into(),
        }
    }
}

impl From<(isize, isize)> for Position {
    fn from((x, y): (isize, isize)) -> Self {
        Position { x, y }
    }
}

impl Position {
    fn sqr_len(&self) -> isize {
        self.x * self.x + self.y + self.y
    }
}

impl FPosition {
    fn sqr_len(&self) -> SmallFraction {
        self.x * self.x + self.y + self.y
    }
}

impl Mul<FPosition> for SmallFraction {
    type Output = FPosition;

    fn mul(self, rhs: FPosition) -> Self::Output {
        FPosition {
            x: self * rhs.x,
            y: self * rhs.y,
        }
    }
}

impl Sum for Position {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        iter.fold(Position { x: 0, y: 0 }, |acc, pos| acc + pos)
    }
}

impl Add for Position {
    type Output = Position;

    fn add(self, rhs: Self) -> Self::Output {
        Position {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

#[derive(Debug)]
struct Lattice {
    sites_x: usize,
    sites_y: usize,
    start_x: usize,
    start_y: usize,
}

impl Lattice {
    fn new(x_size: usize, y_size: usize) -> Self {
        Lattice {
            sites_x: x_size - 1,
            sites_y: y_size - 1,
            start_x: x_size / 2,
            start_y: y_size / 2,
        }
    }

    fn starting_point(&self) -> Position {
        (self.start_x as isize, self.start_y as isize).into()
    }

    fn get_size(&self) -> (usize, usize) {
        (self.sites_x, self.sites_y)
    }

    fn adjacent(&self, &Position { x, y }: &Position) -> Vec<Position> {
        let mut avec = Vec::with_capacity(4);
        if x < self.sites_x as isize {
            // not on the right edge
            avec.push(Position { x: x + 1, y });
        };
        if x > 0 {
            // not on the left edge
            avec.push(Position { x: x - 1, y });
        };
        if y < self.sites_y as isize {
            // not on the bottom edge
            avec.push(Position { x, y: y + 1 });
        };
        if y > 0 {
            // not on the top edge
            avec.push(Position { x, y: y - 1 });
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
        }
    }
}

impl Sub for Position {
    type Output = Position;

    fn sub(self, rhs: Self) -> Self::Output {
        Position {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

#[derive(Debug)]
struct Polymer<'a> {
    positions: Vec<Position>,
    weights: Vec<u8>,
    head: Position,
    lattice: &'a Lattice,
    len: usize,
}

impl Polymer<'_> {
    fn weight(&self, length: usize) -> BigUint {
        self.weights
            .iter()
            .take(length)
            .fold(BigUint::new(vec![1]), |acc, new| acc * new)
    }
}

enum GrowthResult {
    Success,
    Cornered,
}

trait Observable {
    fn value(polymer: &Polymer) -> Fraction;

    fn weighted(polymers: &[Polymer], length: usize) -> Fraction {
        let numerator: Fraction = polymers
            .par_iter()
            .filter(|poly| poly.len >= length)
            .map(|poly| Self::value(poly) * poly.weight(length).clone())
            .sum();
        let denominator: BigUint = polymers
            .par_iter()
            .filter(|poly| poly.len >= length)
            .map(|poly| poly.weight(length).clone())
            .sum();
        numerator / denominator
    }
}

struct Gyration;
impl Gyration {
    fn pdif_center(
        &Position { x: px, y: py }: &Position,
        &FPosition { x: cx, y: cy }: &FPosition,
    ) -> Fraction {
        let diff_x = Into::<SmallFraction>::into(px) - cx;
        let diff_y = Into::<SmallFraction>::into(py) - cy;
        let sdiff = diff_x * diff_x + diff_y * diff_y;

        sdiff.into_fraction::<BigInt>()
    }
}

impl Observable for Gyration {
    fn value(polymer: &Polymer) -> Fraction {
        let pref =
            Into::<SmallFraction>::into(1) / Into::<SmallFraction>::into(polymer.positions.len());
        let center: FPosition =
            pref * FPosition::from(polymer.positions.iter().map(|&pos| pos).sum::<Position>());
        pref.into_fraction::<BigInt>()
            * polymer
                .positions
                .iter()
                .map(|particle| Self::pdif_center(particle, &center))
                .sum::<Fraction>()
    }
}

struct EndToEnd;
impl Observable for EndToEnd {
    fn value(polymer: &Polymer) -> Fraction {
        (polymer.head - polymer.positions[0]).sqr_len().into()
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
        self.weights.push(valid_adj.len() as u8);
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
            weights: Vec::new(),
            lattice,
            len: 0,
        }
    }

    fn grow_till_end(&mut self) -> usize {
        let mut step = 0;
        while let GrowthResult::Success = self.grow() {
            step += 1;
        }
        step
    }
}

const WALK_AMOUNT: usize = 10000;
const LATTICE_SIZE: usize = 350;

fn plot_observables(
    e2es: Vec<Fraction>,
    gyros: Vec<Fraction>,
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
                .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
                .expect("there should be a gyro"),
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
    root.present()?;

    Ok(())
}

fn plot_walks(size: (usize, usize), polies: &[Polymer]) -> Result<(), Error> {
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
        .build_cartesian_2d(0..size.0, 0..size.1)?;

    chart.configure_mesh().draw()?;

    // draw each polymer
    for walk in (0..WALK_AMOUNT).progress() {
        let part = walk as f32 / WALK_AMOUNT as f32;
        let color = ViridisRGBA::get_color(part);
        if WALK_AMOUNT < 100 {
            println!(
                "walk {walk} has a weights of {:?} and achieved length {}",
                polies[walk].weights, polies[walk].len
            );
        }
        chart.draw_series(LineSeries::new(
            polies[walk]
                .positions
                .iter()
                .map(|&Position { x, y }| (x as usize, y as usize)),
            color.stroke_width(3),
        ))?;
    }
    root.present()?;

    Ok(())
}

fn main() -> Result<(), Error> {
    // initialize the lattice
    let lattice = Lattice::new(LATTICE_SIZE, LATTICE_SIZE);

    // create and walk all of the polymers
    let mut polies: Vec<Polymer> = (0..WALK_AMOUNT).map(|_| Polymer::new(&lattice)).collect();
    let bar = ProgressBar::new(WALK_AMOUNT as u64);
    for poly in &mut polies {
        poly.grow_till_end();
        bar.inc(1);
    }
    bar.finish_with_message("finished generating random walks");

    plot_walks(lattice.get_size(), &polies)?;

    let max_poly_len = polies
        .par_iter()
        .max_by(|a, b| a.len.cmp(&b.len))
        .expect("there is a polymer")
        .len
        - 1;
    let e2es: Vec<_> = (0..max_poly_len)
        .progress()
        .map(|len| EndToEnd::weighted(&polies, len))
        .collect();
    let gyros: Vec<_> = (0..max_poly_len)
        .progress()
        .map(|len| Gyration::weighted(&polies, len))
        .collect();

    for (idx, (gyro, e2e)) in gyros.iter().zip(e2es.iter()).enumerate() {
        println!("{} gyro: {:.2}, e2e: {:.2}", idx, gyro, e2e);
    }

    plot_observables(e2es, gyros, max_poly_len)?;

    Ok(())
}
