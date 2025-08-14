#include <iostream>
#include <string> // string
#include <sys/stat.h>  // mkdir
#include <unistd.h>     // chdir

#include "athena.hpp"

//----------------------------------------------------------------------------------------
//! \fn void ChangeRunDir(const char *pdir)
//  \brief change to input run directory; create if it does not exist yet

void ChangeRunDir(const std::string dir) {
  if (dir.empty()) return;

  mkdir(dir.c_str(), 0775);
  if (chdir(dir.c_str())) {
    std::cout << "### FATAL ERROR in " << __FILE__ << " at line " << __LINE__ << std::endl
              << "Cannot cd to directory '" << dir << "'";
    exit(EXIT_FAILURE);
  }

  return;
}
