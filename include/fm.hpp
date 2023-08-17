#pragma once

#include <cinttypes>

namespace fm {
  enum class type { vec, iso, mat, sym, diag, skew, ortho };
}

#include "types/vec.hpp"
#include "types/iso.hpp"
#include "types/mat.hpp"
#include "types/sym.hpp"
#include "types/diag.hpp"
#include "types/skew.hpp"
#include "types/ortho.hpp"

#include "operations/dot.hpp"
#include "operations/outer.hpp"
#include "operations/diag_iso.hpp"
#include "operations/linear_solve.hpp"