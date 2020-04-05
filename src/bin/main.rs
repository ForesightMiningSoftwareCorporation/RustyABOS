use abos::{ABOSGrid};
extern crate rand;
//use rand::prelude::*;

fn main() {
    let mut points: Vec<Vec<f64>> = vec![];
    for ii in 0..3 {
        //making f64 then converting seemse better than casting f64 3 times
        let iif = 1.0 + ii as f64;
        let new_point:Vec<f64> = vec!(iif , iif*2.0, iif*3.0);
        //let new_point:Vec<f64> = vec!(ii , (ii*2) as f64, (ii*3) as f64);
        points.push(new_point);
    }
    
    let mut test = ABOSGrid::new(points,30.0, 0);

    //iteration cycle
    // 1. Filtering points XYZ, specification of the grid, computation of the matrices NB and
    // K, Z→DZ, 0→DP
    //test.init_distance_point_matrixes(); //prepare nb/k/kmax //jmw is this done every loop? if not I've moved it to new function
    // 2. Per partes constant interpolation of values DZ into the matrix P
    test.per_parts_constant_interpolation();
    // 3. Tensioning and smoothing of the matrix P
    ABOSGrid::tension_loop(test.n, test.i1 as usize, test.j1 as usize, &mut test.p, & test.k);
    // 4. P+DP→P
    // 5. ( , )
    // i Xi Yi Z − f →DZi
    // 6. If the maximal difference max { DZi
    //     ,i=1,, n } does not exceed defined precision, the algorithm is finished
    // 7. P→DP, continue from step 2 again (= start the next iteration cycle)
    test.output_all_matrixes();
}
