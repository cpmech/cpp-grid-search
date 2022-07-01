const NNODE: usize = 3;

fn print_x<'a, F>(num_triangle: usize, get_x: F)
where
    F: Fn(usize, usize) -> &'a [f64],
{
    for t in 0..num_triangle {
        for m in 0..NNODE {
            let x = get_x(t, m);
            println!("triangle # {}: x = {:?}", t, x);
        }
    }
}

fn main() {
    // [num_triangle][nnode=3][ndim=2]
    let triangles = vec![
        vec![vec![0.0, 0.0], vec![1.0, 0.0], vec![0.0, 1.0]],
        vec![vec![1.0, 0.0], vec![1.2, 1.5], vec![0.0, 1.0]],
    ];

    // closure that returns the coordinates of cell's point i
    let get_x = |t: usize, m: usize| &triangles[t][m][..];

    // print data
    print_x(triangles.len(), get_x);
}
