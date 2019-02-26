


__kernel void pl_load_grid(
                           __global double2 *grid,
                           double2 origin,
                           double2 size
) {
    int i = get_global_id(0);
    int j = get_global_id(1);
    
    int i_size = get_global_size(0);
    int j_size = get_global_size(1);
    
    grid[i*j_size+j].x = origin.x + size.x * j / (j_size - 1);
    grid[i*j_size+j].y = origin.y + size.y * i / (i_size - 1);
}

__kernel void pl_cartesian_apply_affine_transform_2d(
    __global double16 *xy_in,
    __global double16 *xy_out,
    double8 matrix
) {
    int i = get_global_id(0);

    double8 x = xy_in[i].even;
    double8 y = xy_in[i].odd;
    
    xy_out[i].even = matrix.s0 * x + matrix.s1 * y + matrix.s2;
    xy_out[i].odd = matrix.s3 * x + matrix.s4 * y + matrix.s5;
}
