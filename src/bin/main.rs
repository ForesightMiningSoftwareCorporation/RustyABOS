use abos::{ABOSGrid};
extern crate rand;
use rand::prelude::*;

fn main() {
    let mut points: Vec<Vec<f64>> = vec![];
    let mut rng = rand::thread_rng();
    for _ in 0..10 {
        let new_point:Vec<f64> = vec!(rng.gen_range(200.0,340.0),rng.gen_range(100.0,200.0),rng.gen_range(0.0,10.0));
        points.push(new_point);
    }

    let _test = ABOSGrid::new(points,10000.0, 0);
}
//
// The size of the regular rectangular grid is set according to the following points:
// 1. The greater side of the rectangular domain D is selected, i.e. greater number of
// x21= x2−x1 and y21= y2−y1 . Without loss of generality we can assume
// that x21 is greater.
// 2. The minimal grid size is computed as i0=round  x21/ Dmc  , where Dmc is the
// minimal Chebyshev distance between pairs of points XYZ:
// Dmc=min {max {∣Xi−X j
// ∣,∣Yi−Y j
// ∣} ; i≠ j ∧ i , j=1,, n}
// 3. The optimal grid size is set as: i1=max { k⋅i0 ; k=1, ,5 ∧ k⋅i0Filter }
// 4. The second size of the grid is: j1=round  y21/ x21⋅i1−11
// The presented procedure ensures that the difference between Dx and Dy is minimal i.e. the
// regular rectangular grid is as close to a square grid as possible.

