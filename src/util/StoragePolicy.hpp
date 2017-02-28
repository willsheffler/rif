#ifndef INCLUDED_scheme_util_StoragePolicy_HH
#define INCLUDED_scheme_util_StoragePolicy_HH

#include "types.hpp"

namespace scheme {
namespace util {

///////////////////////////////////////////////////////////////////////////////////////////////
/// storage policies
///////////////////////////////////////////////////////////////////////////////////////////////

///@brief Storage Policy Class, store by value
///@details this is a sortof DUMMY baseclass to define the StoragePolicy concept
///  also a placeholder for doxygen
///@tparam Value type of value that is stored
template <class Value>
class StoragePolicy {
  // StoragePolicy(Value const &);
  ///@return reference to the stored value
  // Value const & value() const;
 protected:
  ///@brief set the value
  void set_stored_value(Value const &v);
  ///@return nonconst reference to stored Value
  // Value & nonconst_value();
};

/// @brief Storage Policy Class, store by value
/// @tparam Value ValueType stored
template <class Value>
struct StoreValue {
  StoreValue() : value_() {}
  StoreValue(Value const &v) : value_(v) {}
  ///@return const reference to stored Value
  Value const &value() const { return value_; }

 protected:
  ///@brief set the value
  void set_stored_value(Value const &v) { value_ = v; }
  ///@return nonconst reference to stored Value
  Value &nonconst_value() { return value_; }
  Value value_;
  ~StoreValue() {}
};

/// @brief Storage Policy Class, store nothing
template <class Value>
struct StoreNothing {
  StoreNothing() {}
  // StoreNothing(Value const & v) : value_(v){}
  ///@return const reference to stored Value
  // Value const & value() const { return value_; }
 protected:
  ///@brief set the value
  void set_stored_value(Value const &v) {}
  ///@return nonconst reference to stored Value
  // Value & nonconst_value() { return value_; }
  ~StoreNothing() {}
};

/// @brief Store-by-pointer policy
/// @tparam Value ValueType stored
/// @note addes the ability to set the pointer
template <class Value>
struct StorePointer {
  ///@return const reference to stored Value
  Value const &value() const { return *value_; }
  /// @brief switch the pointer the this policy manages
  /// @param new_pointer
  void set_pointer(Value *new_pointer) { value_ = new_pointer; }

 protected:
  ///@brief set the value
  void set_stored_value(Value const &v) { *value_ = v; }
  ///@return nonconst reference to stored Value
  Value &nonconst_value() { return *value_; }
  Value *value_;
  ~StorePointer() {}
};

/// @brief Store-by-pointer policy
/// @tparam Value ValueType stored
/// @note addes the ability to set the pointer
template <class Value>
struct StoreSharedPointer {
  ///@return const reference to stored Value
  Value const &value() const { return *value_; }
  /// @brief switch the pointer the this policy manages
  /// @param new_pointer
  void set_pointer(Value *new_pointer) { value_ = new_pointer; }

 protected:
  ///@brief set the value
  void set_stored_value(Value const &v) { *value_ = v; }
  ///@return nonconst reference to stored Value
  Value &nonconst_value() { return *value_; }
  shared_ptr<Value> value_;
  ~StoreSharedPointer() {}
};
}
}

#endif
