NAME          moore90_2
ROWS
 L  R0001
 L  R0002
 L  R0003
 N  R0004
COLUMNS
    INT1      'MARKER'                 'INTORG'
    C0001     R0001     -1
    C0001     R0002     -1
    C0001     R0003     2.5
    C0001     R0004     1
    C0002     R0001     2.5
    C0002     R0002     -2.5
    C0002     R0003     1
    C0002     R0004     2
    INT1END   'MARKER'                 'INTEND'
RHS
    B         R0001     3.75
    B         R0002     -3.75
    B         R0003     8.75
BOUNDS
 UP BOUND     C0001     3
 LO BOUND     C0002     1
 UP BOUND     C0002     2
ENDATA
