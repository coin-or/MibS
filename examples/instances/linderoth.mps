NAME          linderoth   
ROWS
 N  OBJROW
 L  R0000000
 L  R0000001
 L  R0000002
 L  R0000003
 G  R0000004
COLUMNS
    C0000000  OBJROW     -4.           R0000000  1.          
    C0000000  R0000002  -4.             R0000003  1.          
    C0000001  OBJROW    1.             R0000000   -1.        
    C0000001  R0000002   1.           R0000004  1.          
    C0000002  OBJROW     -2.           R0000001  1.          
    C0000002  R0000003  1.          
    C0000003  OBJROW    2.             R0000001  4.          
    C0000003  R0000004  1.             
    C0000004  OBJROW     -6.           R0000000  1.          
    C0000004  R0000001   -3.        
    C0000005  OBJROW    3.             R0000000   -4.        
    C0000005  R0000001  1.             R0000002  1.          
RHS
    RHS       R0000000  7.             R0000001  4.          
    RHS       R0000002  2.             R0000003  1.          
    RHS       R0000004  1.            
BOUNDS
 BV BOUND     C0000000  1.          
 BV BOUND     C0000001  1.          
 BV BOUND     C0000002  1.          
 BV BOUND     C0000003  1.                             
 UI BOUND     C0000004  1e+30
 UI BOUND     C0000005  1e+30
ENDATA
