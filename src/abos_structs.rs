//! Input, Output, and Intermmediate Structs for ABOS
#![warn(missing_docs)]
extern crate nalgebra as na;

use nalgebra::{DVector, Dim, Dynamic, MatrixMN, U3};

pub const INFINITY: f64 = 1.0f64 / 0.0f64;

///Input points and Algorithm inputs
///Grid specified by grid files, or autoset by function call 
pub struct ABOSInputs {
    //User Inputs
    pub linear_tensioning_degree: i8,
    pub filter: f64,
    //INPUT resolution parameter
    pub points: Vec<Vec<f64>>,
    pub q_smooth: f64,
    //user input degree of smoothing default 0.5
    pub grid_enlargement: i32
}

///Outputs to describe the calculated grid
pub struct ABOSOutputs {
    //TODO, should we add in the filtered xyz points computed on?
    /// 2d Matrix of Z values at cell bottom left corner
    pub p: MatrixMN<f64, Dynamic, Dynamic>,
    pub x_min: f64,
    pub y_min: f64,
    pub dx: f64,
    pub dy: f64
}

///Variables of the ABOS algorithm
pub struct ABOSMutable {
    //TODO add p_prior for parallelization

    /// Elements of the elevation matrix [i,j] size
    pub p: MatrixMN<f64, Dynamic, Dynamic>,
    /// Intermediate matrix which holds previous P value or zero
    pub(crate) dp: MatrixMN<f64, Dynamic, Dynamic>,
    /// Matrix of weights used in smoothing, same sze as p & dp
    pub(crate) t_smooth: MatrixMN<f64, Dynamic, Dynamic>,
    /// Vector if z coordinates XYZ
    pub dz: DVector<f64>,
}

///Constants of the ABOS algorithm
pub struct ABOSImmutable {
    /// Degree of Linear Tensioning specifed in Art of Surface Interpolation 3.4.2
    pub(crate) degree: i8,
    /// parameter specifed by degree in Art of Surface Interpolation 3.4.2
    pub(crate) r: usize,
    /// parameter specifed by degree in Art of Surface Interpolation 3.4.2
    pub(crate) l: f64,
    /// Points are potentially filtered to prevent multiple points in the same cell
    pub(crate) xyz_points: MatrixMN<f64, Dynamic, U3>,
    //// Min x
    pub(crate) x1: f64,
    //// Max x
    pub(crate) x2: f64,
    /// Min y
    pub(crate) y1: f64,
    /// Max y
    pub(crate) y2: f64,
    /// Min z
    pub(crate) _z1: f64,
    /// Max z
    pub(crate) _z2: f64,
    
    //ABOS Calculated Parameters
    /// xsize of grid
    pub(crate) i1: i32,
    /// Ysize of grid
    pub(crate) j1: i32,
    /// Cell size in x
    pub(crate) dx: f64,
    /// Cell size in y
    pub(crate) dy: f64,
    /// Matrix of nearest points on grid. Containing indexes to nearest point in the XYZ array
    pub nb: MatrixMN<usize, Dynamic, Dynamic>,
    /// Vector if z coordinates XYZ
    pub(crate) z: DVector<f64>,
    
    /// Grid distance of each grid to the point indexed in NB
    /// k is max, second u is x dist, v is y dist
    pub k_u_v: MatrixMN<(usize, usize, usize), Dynamic, Dynamic>,

    /// Maximal element of matrix K
    pub(crate) k_max: usize,
    /// Whether the xy points were swapped in pre process so swapped post process
    pub(crate) xy_swaped: bool,
    /// User input smoothing parameter
    pub(crate) q_smooth: f64,
}

impl ABOSImmutable {
    ///convert indexes to xy points.  Point at cell bottom left point.
    pub(crate) fn indexes_to_position(&self, row_index: &usize, col_index: &usize) -> [f64; 2] {
        let x_coordinate = self.x1 + self.dx * (*row_index as f64);
        let y_coordinate = self.y1 + self.dy * (*col_index as f64);

        [x_coordinate, y_coordinate]
    }
}
