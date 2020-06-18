//
// Created by Moritz Gnisia on 05.04.20.
//

#ifndef CFDLAB_ENUMS_H
#define CFDLAB_ENUMS_H

enum class velocity_type {
    U,
    V
};

enum class border_position {
    TOP,
    BOTTOM,
    LEFT,
    RIGHT
};

enum class matrix_selection {
    ROW,
    COLUMN
};

enum class read_type {
    INT,
    DOUBLE
};

enum class ID {
    A,
    B,
    C,
    D
};

#endif //CFDLAB_ENUMS_H
