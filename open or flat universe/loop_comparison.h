#ifndef LOOP_COMPARISON_H
#define LOOP_COMPARISON_H

//Function objects to serve as loop conditions for particle propagation
struct Lattice::compare_block_t
{
    const Particle &target_p;
    compare_block_t(const Particle &p) : target_p(p) {}
    bool operator()(const Particle &p) const {return block::compare_time(target_p.block_p, p.block_p) > 0;}
};
struct Lattice::compare_r
{
    const Particle &target_p;
    compare_r(const Particle &p) : target_p(p) {}
    bool operator()(const Particle &p) const {return block::compare_rad(target_p.block_p, p.block_p) > 0;}
};
struct Lattice::compare_face
{
    const block::Face target_face;
    compare_face(block::Face x) : target_face(x) {}
    bool operator()(const Particle &p) const {return target_face != p.face;}
};
struct Lattice::change_traj_dir
{
    enum sign {neg=-1, pos=1};
    const sign s;
    
    change_traj_dir(sign in_sign) : s(in_sign) {}
    bool operator()(const Particle &p) const {return s*p.traj[0] < 0;}
};

#endif
