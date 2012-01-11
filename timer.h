
#include <unistd.h>
#include <cstdlib>
#include <time.h>
#include <list>

class split {

public:
  std::string label;
  time_t when;

  split(std::string l);
};

class timer {
private:
  std::string title;
  std::list<split> splits;
  float timesub(time_t t1, time_t t2);
public:
  timer(std::string title);
  void addsplit(std::string l);
  void print();
};
