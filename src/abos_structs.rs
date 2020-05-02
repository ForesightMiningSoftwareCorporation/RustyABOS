extern crate nalgebra as na;

use nalgebra::{DVector, Dynamic, MatrixMN, U3};

pub const INFINITY: f64 = 1.0f64 / 0.0f64;

pub struct ABOSInputs {
    //User Inputs
    pub degree: i8,
    pub filter: f64,
    //INPUT resolution parameter
    pub points: Vec<Vec<f64>>,
    pub q_smooth: f64,
    //user input degree of smoothing default 0.5
}

pub struct ABOSMutable {
    //size of grid on y
    pub p: MatrixMN<f64, Dynamic, Dynamic>,
    //elements of the elevation matrix [i,j] size
    pub(crate) dp: MatrixMN<f64, Dynamic, Dynamic>,
    //matrix of weights used in smoothing, same sze as p & dp
    pub(crate) t_smooth: MatrixMN<f64, Dynamic, Dynamic>,
    //vector if z coordinates XYZ
    pub dz: DVector<f64>,
}

pub struct ABOSImmutable {
    pub(crate) degree: i8,
    pub(crate) r: usize,
    pub(crate) l: f64,
    pub(crate) xyz_points: MatrixMN<f64, Dynamic, U3>,
    //Derived Points
    pub(crate) x1: f64,
    //minx
    pub(crate) x2: f64,
    //maxx
    pub(crate) y1: f64,
    //miny
    pub(crate) y2: f64,
    //maxy
    pub(crate) _z1: f64,
    //min z
    pub(crate) _z2: f64,
    //max z
    //ABOS Calculated Parameters
    pub(crate) dmc: f64,
    /// TODO WHAT IS THIS FOR??
    //minimal chebyshev distance
    pub(crate) i1: i32,
    //xsize of grid
    pub(crate) j1: i32,
    //ysize of grid
    pub(crate) dx: f64,
    //size of grid on x
    pub(crate) dy: f64,
    //auxiliary matrix same size as P
    pub nb: MatrixMN<usize, Dynamic, Dynamic>,
    //[i x j]
    pub(crate) z: DVector<f64>,
    //auxiliary vector same size as Z
    pub k_u_v: MatrixMN<(usize, usize, usize), Dynamic, Dynamic>,
    //first one is max, second one is x dist, third one is y dist
    // Grid distance of each grid to the point indexed in NB
    pub(crate) k_max: usize,
    //maximal element of matrix K
    pub(crate) _rs: f64,
    // Resolution of map
    pub(crate) _xy_swaped: bool,
    //whether the xy points were swapped
    pub(crate) q_smooth: f64,
    //user input smoothing parameter
}

impl ABOSImmutable {
    pub(crate) fn indexes_to_position(&self, row_index: &usize, col_index: &usize) -> [f64; 2] {
        let x_coordinate = self.x1 + self.dx * (*row_index as f64);
        let y_coordinate = self.y1 + self.dy * (*col_index as f64);

        [x_coordinate, y_coordinate]
    }
}
