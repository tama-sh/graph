#include "graph.h"

std::ostream& operator<<(std::ostream &os, const Edge &e)
{
    os << e.source << " " << e.target;
    return os;
}