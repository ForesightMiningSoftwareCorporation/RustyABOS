use abos::{abos_run_auto_grid, abos_run_manual_grid};
use abos::abos_structs::{ABOSAutoGridInputs, ABOSManualGridInputs, ABOSOutputs};
use abos::io_system::{import_points_csv, output_grd_file};

extern crate rand;
//use rand::prelude::*;

fn main() {
    let mut points: Vec<Vec<f64>> = vec![];

    let ok = import_points_csv(&mut points, "./testFiles/inputs/xyz.csv");
    if ok.is_err(){ return}
    // for ii in 0..3 {
    //     //making f64 then converting seemse better than casting f64 3 times
    //     let iif = 1.0 + ii as f64;
    //     let new_point: Vec<f64> = vec![iif, iif * 2.0, iif * 3.0];
    //     //let new_point:Vec<f64> = vec!(ii , (ii*2) as f64, (ii*3) as f64);
    //     points.push(new_point);
    // }

    let inputs = ABOSAutoGridInputs {
        linear_tensioning_degree: 0,
        filter: 100.0,
        points,
        q_smooth: 0.5,
        grid_enlargement: 5
    };

    let outputs: ABOSOutputs = abos_run_auto_grid(&inputs);
    output_grd_file(outputs, "./testFiles/inputs/out_auto.grd");

    let mut points: Vec<Vec<f64>> = vec![];

    let ok = import_points_csv(&mut points, "./testFiles/inputs/xyz.csv");
    if ok.is_err(){ return}
    let inputs = ABOSManualGridInputs {
        linear_tensioning_degree: 0,
        filter: 100.0,
        points,
        q_smooth: 0.5,
        grid_enlargement: 5,
        x_min: 100.0,
        dx: 100.0,
        nx: 10,
        y_min: 100.0,
        dy: 100.0,
        ny: 10
    };
    let outputs = abos_run_manual_grid(&inputs).expect("Initialization failed, invalid inputs");
    output_grd_file(outputs, "./testFiles/inputs/out_manual.grd");
    
}
