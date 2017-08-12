#ifndef __AXIOM_HPP__
#define __AXIOM_HPP__

#include "common.hpp"

namespace axiom {

class Axiom {
 public:
  ~Axiom();
  static Axiom& Get();
  enum Derivator { CPU, GPU };

#ifndef CPU_ONLY
  inline static cublasHandle_t cublas_handle() { return Get().cublas_handle_; } 
#endif

  inline static Derivator mode()                { return Get().mode_; }
  inline static void set_mode(Derivator mode)   { Get().mode_ = mode; }
  inline static bool multiprocess()             { return Get().multiprocess_; }
  inline static void set_multiprocess(bool val) { Get().multiprocess_ = val; }

  static void SetDevice(const int device_id);
  // Prints the current GPU status.
  static void DeviceQuery();
  // Check if specified device is available
  static bool CheckDevice(const int device_id);
  // Search from start_id to the highest possible device ordinal,
  // return the ordinal of the first available device.
  static int FindDevice(const int start_id = 0);

 protected:
#ifndef CPU_ONLY
  cublasHandle_t cublas_handle_;
#endif
  Derivator mode_;
  bool multiprocess_;
 private:
  /* The private constructor to avoid duplicate instantiation */
  Axiom();
  DISABLE_COPY_AND_ASSIGN(Axiom);

}; /* class Axiom */ 
}  /* namespace axiom */
#endif /* __AXIOM_HPP__ */
