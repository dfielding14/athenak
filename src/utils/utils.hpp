#ifndef UTILS_UTILS_HPP_
#define UTILS_UTILS_HPP_

#include <string>

void ShowConfig();
void ChangeRunDir(const std::string dir);
int CreateMPITag(int lid, int buff_id, int phys_id);

#endif // UTILS_UTILS_HPP_
