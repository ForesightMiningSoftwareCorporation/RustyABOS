///Crate to generate a surface with the ABOS algorithm proposed in Art of Surface Interpolation byM gr. Miroslav Dressler, Ph.D.
/// http://m.dressler.sweb.cz/ABOS.htm
/// 

extern crate approx; // For the macro relative_eq! in nalgebraImportCheck
extern crate nalgebra as na;
use na::{DMatrix, DVector};

// pub struct abosParameters {
//     zmin
// }


fn get_chebyshev_distance(x1 : f64, x2 : f64, y1 : f64, y2 : f64) -> f64{
    if (x1 - x2).abs() > (y1 - y2).abs() {
        (x1 - x2).abs()
    }
    else {
        (y1 - y2).abs()
    }
}
//mutating matrix example
// let mut a = Matrix2x3::new(1, 2, 3, 4, 5, 6);
// for (i, mut row) in a.row_iter_mut().enumerate() {
// row *= (i + 1) * 10;
// }
// let expected = Matrix2x3::new(10, 20, 30, 80, 100, 120);


//TODO make closure 
fn get_min_chebyshev_distance<T: na::Scalar>(X : DVector<T>, Y :DVector<T>) -> f64{
    let min_chebyshev = std::f64::MAX;
    //print!("{:?}, {:?}, {:?}", X, Y, min_chebyshev);
    //loop through vectors
    // let pairMin = |i: f64| -> get_chebyshev_distance(x1: f64, x2: f64, y1: f64, y2: f64);
    // min_chebyshev
    0.0
}



pub fn dynamicMatrixCheck(){
    let (mut x1, mut x2, mut y1, mut y2, mut z1, mut z2) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0); //x & y min(1) and max(2)
    let i1 = 4;
    let j1 = 5;
    let numPoints = 10;

    let mut P  = DMatrix::from_element(i1, j1, 0.0);
    let mut DP  = DMatrix::from_element(i1, j1, 0.0);

    let Z = DVector::from_element(numPoints, 0.0);
    let DZ = DVector::from_element(numPoints, 0.0);
    let X = DVector::from_element(numPoints, 0.0);
    let Y = DVector::from_element(numPoints, 0.0);

    let mut K = DMatrix::from_element(i1, j1, 0.0);
    let mut kMax = K.max();
    let mut NB = DMatrix::from_element(i1, j1, 0.0);
    let mut filter = 1.0; // parameter
    let mut RS = get_chebyshev_distance(x1, x2, y1, y2) / filter;

    x1 = X.min();
    x2 = X.max();
    y1 = Y.min();
    y2 = Y.max();
    z1 = Z.min();
    z2 = Z.max();

    let dx = (x2 - x1) / ((i1 - 1) as f64);
    let dy = (y2 - y1) / ((j1 - 1) as f64);


    println!("RS {:?}", RS);
    println!("x1 {:?}", x1);
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }

    #[test]
    fn nalgebraImportCheck() {
        let axis  = na::Vector3::x_axis();
        let angle = 1.57;
        let b     = na::Rotation3::from_axis_angle(&axis, angle);
    
        approx::relative_eq!(b.axis().unwrap(), axis);
        approx::relative_eq!(b.angle(), angle);
    }
}


pub fn hello1() -> () {
    println!("Hello1");
}

pub fn hello2() -> () {
    println!("Hello2");
}

