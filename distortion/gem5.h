#ifndef TAU_GEM5_H
#define TAU_GEM5_H

#if defined(__arm__) || defined(__AARCH64EL__)
extern "C" {
#include <m5op.h>
}
#endif

namespace gem5 {

  inline void ResetStats()
  {
#if defined(__arm__) || defined(__AARCH64EL__)
    m5_reset_stats(0, 0);
#endif
  }
  inline void DumpStats()
  {
#if defined(__arm__) || defined(__AARCH64EL__)
    m5_dump_stats(0, 0);
#endif
  }

  class AutoStats
  {
    AutoStats() {
      ResetStats();
    }

    ~AutoStats() {
      DumpStats();
    }
  };

}

#endif // TAU_GEM5_H

