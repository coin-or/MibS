NAME          moore90
ROWS
 L  R0001
 L  R0002
 L  R0003
 L  R0004
 N  R0005
COLUMNS
    INT1      'MARKER'                 'INTORG'
    C0001     R0001     -25
    C0001     R0002     1
    C0001     R0003     2
    C0001     R0004     -2
    C0001     R0005     -1
    C0002     R0001     20
    C0002     R0002     2
    C0002     R0003     -1
    C0002     R0004     -10
    C0002     R0005     -10
    INT1END   'MARKER'                 'INTEND'
RHS
    B         R0001     30
    B         R0002     10
    B         R0003     15
    B         R0004     -15
BOUNDS
 UP BOUND     C0001     10
 UP BOUND     C0002     5
ENDATA
