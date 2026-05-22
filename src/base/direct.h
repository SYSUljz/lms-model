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

constexpr double GetDirectFactor(Direct8 direct)
{
  switch (direct)
  {
  case Direct8::DownRight:
  case Direct8::DownLeft:
  case Direct8::UpLeft:
  case Direct8::UpRight:
    return 1.4142135623730951;
  default:
    return 1.0;
  }
}