// this is a factory class to build model

#include "./model.cpp"
class ModelBuilder
{
  model _model;

public:
  ModelBuilder &shape(int heigh, int weigh) {};
  ModelBuilder &file_path() {};
  model &build()
  {
    return _model;
  }
};