
#include <unistd.h>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <time.h>
#include <list>

#include "timer.h"

split::split(std::string l) {
  label=l;
  time(&when);
}

float timer::timesub(time_t t1, time_t t2) {
  float sec;
  if (t2 > t1) {
    sec = (float)(t2-t1);
  } else {
    sec = (float)(t1-t2);
  }
  return (float) sec;
}

timer::timer(std::string t) {
  title=t;
  timer::addsplit("initialized");
}

void timer::addsplit(std::string l) {
  split *s=new split(l);
  splits.push_back(*s);
}
  
void timer::print() {
  std::list<split>::iterator s;
  float delta, fromstart;
  s=splits.begin();
  split &first = *s;
  split *last = &(*s);
  s++;
  std::cout << title.c_str() << "\n";
  std::cout << std::setw(10) << "split" << std::setw(10) << "last" << std::setw(10) << "total\n";
  for (; s!=splits.end(); ++s) {
    fromstart = timesub(first.when, s->when);
    delta = timesub(last->when, s->when);
    last=&(*s);
    std::cout << std::setw(10) << s->label.c_str() << "  " <<
      std::setw(10) << std::setprecision(6) << delta <<
      std::setw(10) << std::setprecision(6) << fromstart <<
      "\n";
  }
}

