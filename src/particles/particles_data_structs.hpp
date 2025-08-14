#ifndef PARTICLES_DATA_STRUCTS_HPP_
#define PARTICLES_DATA_STRUCTS_HPP_

//----------------------------------------------------------------------------------------
//! \struct ParticleLocationData
//! \brief data describing location of data for particles communicated with MPI

struct ParticleLocationData {
  int prtcl_indx;   // index in particle array
  int dest_gid;     // GID of target MeshBlock
  int dest_rank;    // rank of target MeshBlock
};

// Custom operators to sort ParticleLocationData array by dest_rank or prtcl_indx
struct {
  bool operator()(ParticleLocationData a, ParticleLocationData b)
    const { return a.dest_rank < b.dest_rank; }
} SortByRank;
struct {
  bool operator()(ParticleLocationData a, ParticleLocationData b)
    const { return a.prtcl_indx < b.prtcl_indx; }
} SortByIndex;

//----------------------------------------------------------------------------------------
//! \struct ParticleMessageData
//! \brief Data describing MPI messages containing particles

struct ParticleMessageData {
  int sendrank;  // rank of sender
  int recvrank;  // rank of receiver
  int nprtcls;   // number of particles in message
  ParticleMessageData(int a, int b, int c) :
    sendrank(a), recvrank(b), nprtcls(c) {}
};

#endif // PARTICLES_DATA_STRUCTS_HPP_

