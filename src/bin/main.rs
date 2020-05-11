use abos::abos_run_autogrid;
use abos::abos_structs::{ABOSInputs, ABOSOutputs};

extern crate rand;
//use rand::prelude::*;

fn main() {
    let mut points: Vec<Vec<f64>> = vec![];
    for ii in 0..3 {
        //making f64 then converting seemse better than casting f64 3 times
        let iif = 1.0 + ii as f64;
        let new_point: Vec<f64> = vec![iif, iif * 2.0, iif * 3.0];
        //let new_point:Vec<f64> = vec!(ii , (ii*2) as f64, (ii*3) as f64);
        points.push(new_point);
    }

    let inputs = ABOSInputs {
        linear_tensioning_degree: 0,
        filter: 100.0,
        points,
        q_smooth: 0.5,
        grid_enlargement: 5
    };

    //println!("inputs.points {:?}", inputs.points);

    let _outputs = abos_run_autogrid(&inputs);
    // test.per_parts_constant_interpolation();
    // ABOSGrid::calculation_loop(&mut test);
}
