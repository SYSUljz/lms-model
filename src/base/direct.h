#pragma once
enum class Direct8
{
  Right,
  DownRight,
  Down,
  DownLeft,
  Left,
  UpLeft,
  Up,
  UpRight
};

template <typename T>
constexpr T GetDirectFactor(Direct8 direct)
{
  switch (direct)
  {
  case Direct8::DownRight:
  case Direct8::DownLeft:
  case Direct8::UpLeft:
  case Direct8::UpRight:
    return static_cast<T>(1.4142135623730951);
  default:
    return T{1};
  }
}