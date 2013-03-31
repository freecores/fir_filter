--  SUBSISTEMA ..........:  
--  TARJETA .............:  

--  FECHA DE CREACIÓN ...:  19 Jul 2007
--  AUTOR ...............:  DIEGO PARDO
--  REVISIÓN ............:  
--  FECHA ULTIMA REVISIÓN:  
--  AUTOR REVISIÓN ......:  

--  TITULO "FILTRO FIR DE COEFICIENTES SIMETRICOS";


-------------------------------------------------------------------------------
-- DESCRIPCION:
-------------------------------------------------------------------------------
--          
--  Filtro FIR de orden variable (orden min = 1) y totalmente parametrizable.
--  Dispone de control de saturacion negativa y positiva. (CONSTANT max_coef= 64).

--  La instanciacion y declaracion de este componente estan al final de este
--  fichero.

--  Para la obtencion de los coeficientes se puede usar el fichero de matlab 
--  "filtros_fir.m". En ese fichero solo es necesario modificar "fs", "orden"
--  y "fc", y directamente devuelve los coeficientes. Tambien se puede ver unas
--  ilustraciones de respuesta a diferentes señales de entrada (escalon, delta 
--  y rampa ascendente).
--  Se pueden introducir coef obtenidos por otros metodos, pero siempre deben
--  ser normalizados a 1 (si no lo estuviesen).
--
-- VELOCIDAD
------------------
--  La velocidad viene marcada por el orden del filtro. El retardo max de salida es
--  igual al orden del filtro + 1 ciclo.
--
-- AREA
------------------
--  La resolucion, como parametro de entrada, marca la precision de las operaciones
--  de filtrado (y por tanto el area que ocupa). A menor resolucion menor area.
--
-- ERROR DE SALIDA
------------------
--  El error maximo a la salida del filtro es igual a: 
-- 
--             -------------------------------------------------
--            | num_coef * [2^(-bits_resol)] * [(2^bits_in)-1)] |
--             -------------------------------------------------
-- (num_coef = orden+1)

-------------------------------------------------------------------------------
-- FUNCIONAMIENTO:
-------------------------------------------------------------------------------

-- x[n]     entrada sincrona positiva
-- h[n]     coeficientes                
-- y[n]     salida sincrona positiva    (procesada en <1 ciclo, sin usar pipeline)

------------------------------------------------------------------------------
-- Operaciones de filtrado:
------------------------------------------------------------------------------

--          coeficentes normalizados    A(0,bits_resol) con resolucion = 2^(-bits_resol)
--          valores x[n]                U(bits_in,0) -> A(bits_in,0)
--
--  y[n] = sum{x[n-k] x h[k]}:
--          productos A(bits_in,0) x A(0,bits_resol) = A(bits_in+0+1,bits_resol)
--          num_coef-1 sumas A(bits_in+1,bits_resol) => y[n] : A(bits_in+1+[num_coef]-1,bits_resol) 
--
------------------------------------------------------------------------------
-- Implementacion:
------------------------------------------------------------------------------
--
--  -FIR en forma directa II-
--
--  Se implementan 3 procesos que operan en paralelo de la siguiente forma:
--
--  proceso1():    realiza operaciones de filtrado con la mitad de coef
--  proceso2():    realiza operaciones de filtrado con la mitad de coef
--  proceso_reg(): implementa el registro de desplazamiento para x[n-k]
--
--  Los filtros FIR son simetricos, por lo que se aprovecha para que cada uno
--  de los procesos que realizan operaciones lo hagan de forma simplificada,

--  p.e. para num_coef=5 y coef positivos: 
--      x[n-3]*h(3) + x[n-1]*h[1] = (x[n-3]+x[n-1])*h[3], ya que h[1]=h[3]
--
---------------
-- 
--  x: productos realizados por proceso1() -> mitad coef intermedios y coef central
--  o: productos realizados por proceso2() -> mitad coef intermedios y coef inicial y final (extremos)
--
--   num_coef             bi             prod proceso1()             prod proceso2()
--
--      2?      x x                             0                           1   bo=b1
--      3       o x o                           1                           1   bo=b2
--      4       o x x o                         1   b1=b2                   1   bo=b3       
--      5       o x x x o                       2   b2=b3                   1   bo=b4       
--      6       o x x x x o                     2   b1=b4, b2=b3            1   bo=b5           
--      7       o o x x x o o                   2   b2=b4                   2   bo=b6, b1=b5 
--      8       o o x x x x o o                 2   b2=b5, b3=b4            2   bo=b7, b1=b6
--      9       o o x x x x x o o               3   b2=b6, b3=b5,           2   bo=b8, b1=b7
--      10      o o x x x x x x o o             3   b2=b7, b3=b6, b4=b5     2   bo=b9, b1=b8
--      11      o o o x x x x x o o o           3   b3=b7, b4=b6            3   bo=b10,b1=b9, b2=b8
--      12      o o o x x x x x x o o o         3   b3=b8, b4=b7, b5=b6     3   bo=b11,b1=b10,b2=b9
--      13      o o o x x x x x x x o o o       4   b3=b9, b4=b8, b5=b7     3   bo=b12,b1=b11,b2=b10
--      14      o o o o x x x x x x o o o o     4   b4=b9, b5=b8, b6=b7     4   bo=b13,b1=b12,b2=b11,b3=b10
--      15      o o o o x x x x x x x o o o o   4   b4=b10,b5=b9, b6=b8     4   bo=b14,b1=b13,b2=b12,b3=b11
--      ...

-------------------------------------------------------------------------------
-- Fixed-Point Arithmetic: An Introduction
-- Randy Yates
-- March 3, 2001 11:52
-------------------------------------------------------------------------------

-- A(a,b): 'a' parte entera sin signo y 'b' parte decimal
--  
--          bits                =   (a+b) + 1 (signo)
--          resolucion (pasos)  =   2^(-b)
--          max value           =   2^a - 2^(-b)
--          min value           =   -2^(-a)
--
--          A(x,y) * A(z,n)         =   A(x+z+1,y+n)    con (x+z+1+y+n) +1 bits (signo)
--          A(x,y) + A(x,y)         =   A(x+1,y)        con (x+1+y) +1 bits (signo)
--          A(x,y) + A(x,y) + A(x,y)=   A(x+2,y)        con (x+2+y) +1 bits (signo)
-- 

------------------------------------------------------------------------------          
--                 CONVERSIONES ENTRE TIPOS (CONSTANTES)
------------------------------------------------------------------------------

--
--      REAL -> INTEGER -> STD_LOGIC_VECTOR -> SIGNED:
--
--          SIGNED(CONV_STD_LOGIC_VECTOR(INTEGER(ROUND(valor_real))))
--
--      SIGNED -> STD_LOGIC_VECTOR -> INTEGER -> REAL:
--
--          ROUND(REAL(CONV_INTEGER(CONV_STD_LOGIC_VECTOR(valor_signed))))
--

--=============================================================================
--============================= PACKAGE LOCAL =================================
--=============================================================================

LIBRARY ieee;
USE ieee.std_logic_1164.all;
USE ieee.std_logic_arith.all;
USE ieee.std_logic_signed.all;
USE ieee.math_real.all;     -- interfaz para introducir los coef en formato real


PACKAGE tipos_filtro IS

    CONSTANT max_coef: NATURAL:= 64; 
    TYPE COEFICIENTES IS ARRAY (max_coef-1 DOWNTO 0) OF REAL;
    FUNCTION real_a_signed (CONSTANT c_real: IN REAL; CONSTANT precision: IN NATURAL) RETURN SIGNED;

END tipos_filtro;


PACKAGE body tipos_filtro IS

    FUNCTION real_a_signed (CONSTANT c_real: IN REAL; CONSTANT precision: IN NATURAL) RETURN SIGNED IS
        -- SIGNED = REAL / 2^(-precision) = REAL * 2^precision, con redondeo IEEE
        CONSTANT RESOLUCION : REAL := "**"(2,REAL(INTEGER(precision)));
    BEGIN
        RETURN (SIGNED(CONV_STD_LOGIC_VECTOR(INTEGER(ROUND ("*"(c_real,RESOLUCION))),precision+1)));
    END FUNCTION;
    
END tipos_filtro;

--=============================================================================
--=========================== ENTITY filtro_FIR ===============================
--=============================================================================

USE WORK.tipos_filtro.all;

LIBRARY ieee;
USE ieee.std_logic_1164.all;
USE ieee.std_logic_arith.all;
USE ieee.std_logic_signed.all;
USE ieee.math_real.all;     -- usada como interfaz para introducir los coef en formato real

--LIBRARY altera;
--USE altera.altera_syn_attributes.all; 


ENTITY filtro_FIR IS
    GENERIC
    (
        bits_in     : NATURAL       := 8;       -- = bits de entrada = bits de salida
        num_coef    : NATURAL       := 32;      -- = N (orden) + 1,  min 2 coef (N = 1)
        bits_resol  : NATURAL       := 32;      -- = bits parte decimal de los coef => resol = 2^(-bits_resol)              
        coeficiente : COEFICIENTES  := (        -- coeficientes bi normalizados (bo, ..., bN)
                0.00475562935382438,
                0.00531548289351449,
                0.00697197720225127,
                0.00965735079759454,
                0.0132617084564161,
                0.017637519982969,
                0.0226056603567445,
                0.0279627439350249,
                0.0334894523768257,
                0.038959515246785,
                0.0441489755173365,
                0.0488453605108728,
                0.05285638268919,
                0.0560178139472758,
                0.0582002109141864,
                0.0593142158191887,
                0.0593142158191887,
                0.0582002109141864,
                0.0560178139472758,
                0.05285638268919,
                0.0488453605108728,
                0.0441489755173365,
                0.038959515246785,
                0.0334894523768257,
                0.0279627439350249,
                0.0226056603567445,
                0.017637519982969,
                0.0132617084564161,
                0.00965735079759454,
                0.00697197720225127,
                0.00531548289351449,
                0.00475562935382438,    
                OTHERS=>0.0)
    );
    PORT
    (
        reset           : IN    STD_LOGIC;                             
        reloj           : IN    STD_LOGIC;                              
        enable          : IN    STD_LOGIC;                              
        xn              : IN    STD_LOGIC_VECTOR(bits_in-1 DOWNTO 0);   --  Dato de entrada sínc
        yn              : OUT   STD_LOGIC_VECTOR(bits_in-1 DOWNTO 0)    --  Dato de salida con control de sat neg
--      yn_fraccional   : OUT   SIGNED (bits_in+1+bits_resol+1+(num_coef-2) DOWNTO 0) -- sin control de sat para depuracion
    );
END filtro_FIR;

ARCHITECTURE arch_filtro_FIR OF filtro_FIR IS
    type ENTRADA_xn IS ARRAY (num_coef-2 DOWNTO 0) OF SIGNED(bits_in DOWNTO 0);     -- A(bits_in,0)     : x[n-k]
    type VALORES_h  IS ARRAY (num_coef-1 DOWNTO 0) OF SIGNED(bits_resol DOWNTO 0);  -- A(0,bits_resol)  : h[n-k]

    SIGNAL REGISTRO_xn  : ENTRADA_xn;       
        
    SIGNAL yn_1, yn_2   : SIGNED (bits_in+1+bits_resol+1+(num_coef-2) DOWNTO 0);
    SIGNAL yn_aux       : SIGNED (bits_in+1+bits_resol+1+(num_coef-2) DOWNTO 0);    -- yn_1 + yn_2

                                                                            
BEGIN

    -- =======================================================================================
    --  PROCESO DE LA MITAD DE LAS OPERACIONES DE FILTRADO 
    --  (NO INCLUYE COEF INICIAL NI FINAL; INCLUYE COEF INTERMEDIO SI "num_coef" ES IMPAR)
    -- =======================================================================================

    proceso1:
    PROCESS (reloj, reset)                      
    
        VARIABLE xn_a_pares : SIGNED (bits_in+1 DOWNTO 0);                          -- para sumas de pares de x[n-k] con coef simetricos
        VARIABLE REGISTRO_h : VALORES_h := (OTHERS => (OTHERS => '0'));             -- para guardar los coef (evita el registro)
        VARIABLE producto   : SIGNED (bits_in+1+bits_resol+1 DOWNTO 0);             -- prod de xn_a_pares por cada coef
        VARIABLE yn_aux     : SIGNED (bits_in+1+bits_resol+1+(num_coef-2) DOWNTO 0);-- sumas parciales
        VARIABLE signos_bi  : STD_LOGIC_VECTOR (1 DOWNTO 0);                        -- signos de los coef simetricos
        VARIABLE signos_xn  : STD_LOGIC_VECTOR (1 DOWNTO 0);                        -- Ssignos de los pares de x[n-k]


        -- estas variables sirven para definir cuantas operaciones el proceso para hacer las 
        -- operaciones de filtrado. El proceso comienza en el coef "indice1" y realiza operaciones 
        -- con un numero  de coeficientes dado por "margen". No se usan ctes porque solo se asignan 
        -- si num_coef supera cierto valor (sino darian error por la formula aplicada).
        VARIABLE indice1    : NATURAL;
        VARIABLE margen     : NATURAL;

    BEGIN
        IF (reset = '1') THEN
            yn_1    <= (OTHERS => '0');
                    
        ELSIF (reloj'event AND reloj='1') THEN
            IF (enable = '1') THEN
            
                -- se usa una variable para evitar el registro de h[n], REGISTRO_h(0) = bN 
                FOR ind IN 0 TO (num_coef-1) LOOP
                    REGISTRO_h(ind) := real_a_signed (coeficiente(max_coef-num_coef+ind),bits_resol);
                END LOOP;           
    
                ----------------------------------------------------------------                                                
                -- operaciones filtrado (num_coef marca el num de prod y sumas)
                ----------------------------------------------------------------
                
                -- en este ciclo aun no se ha aplicado x[n-k] al registro de valores de x[n] 
                                
                -- como los filtros FIR son simetricos, se aprovecha esto para realizar
                -- la mitad de multiplicaciones (considerando los signos): primero se suman entre
                -- si los x[n-k] que corresponden a operaciones con coeficientes simetricos y
                -- luego se le aplica el signo (se extiende siempre el signo '0' al registro de
                -- x[n-k] porque x[n] es siempre positivo)
                yn_aux := (OTHERS => '0');
                
                -- operaciones con los coef simetricos (ni con el coef central ni con el inicial o final):
                
                -- (en caso de que el num de coef sea 3 coinciden los coef simetricos intermedios
                -- con el inicial y final, q son caso especial -> debemos evitar confundirlos)
                        
                IF (num_coef > 3) THEN
                
                    indice1 := (num_coef/2)+(num_coef MOD 2);
                    margen  := ((num_coef-1)-(num_coef/2))/2 - (num_coef MOD 2);
                
                    FOR i IN indice1 TO (indice1 + margen) LOOP
                        
                        -- signos de los bi simetricos
                        signos_bi(1):= REGISTRO_h(i)(bits_resol);
                        signos_bi(0):= REGISTRO_h(num_coef-1-i)(bits_resol);
                    
                        -- signos de los x[n-k] con coef simetricos para extension de signo
                        signos_xn(1):= REGISTRO_xn(i)(REGISTRO_xn(i)'left);
                        signos_xn(0):= REGISTRO_xn(num_coef-1-i)(REGISTRO_xn(num_coef-1-i)'left);
                    
                        CASE signos_bi IS
                            WHEN "00" =>    -- FIR tipo I - II
                                xn_a_pares:= "+"(signos_xn(1)&REGISTRO_xn(i), signos_xn(0)&REGISTRO_xn(num_coef-1-i));
                            WHEN "01" =>    -- FIR tipo III - IV
                                xn_a_pares:= "+"(signos_xn(1)&REGISTRO_xn(i), -(signos_xn(0)&REGISTRO_xn(num_coef-1-i)));
                            WHEN "10" =>    -- FIR tipo III - IV
                                xn_a_pares:= "+"(-(signos_xn(1)&REGISTRO_xn(i)), signos_xn(0)&REGISTRO_xn(num_coef-1-i));
                            WHEN OTHERS =>  -- FIR tipo I - II
                                xn_a_pares:= "+"(-(signos_xn(1)&REGISTRO_xn(i)), -(signos_xn(0)&REGISTRO_xn(num_coef-1-i)));
                        END CASE;               
            
                        producto:= "*"(xn_a_pares, ABS(REGISTRO_h(i)));
                        yn_aux := "+"(yn_aux, producto);
                    END LOOP;       
                END IF;
                
                --------------------------------------------------------------- 
                -- Operacion con el coef central (si num_coef es impar)
                ---------------------------------------------------------------
                
                -- El coef central (si num_coef es impar) no es simetrico con 
                -- ningun otro coef
                
                IF (num_coef MOD 2)/=0 THEN
	                signos_xn(1):= REGISTRO_xn(num_coef/2)(REGISTRO_xn(num_coef/2)'left);                
                    producto    := "*"(signos_xn(1)&REGISTRO_xn(num_coef/2), REGISTRO_h(num_coef/2));
                    yn_aux      := "+"(yn_aux, producto);
                END IF; 
                            
                --------------------------------------------------------------- 
                -- SALIDA FRACCIONAL SIN CONTROL DE SATURACION
                ---------------------------------------------------------------

                -- yn_1 agrupa los bits de la parte entera y decimal (bits_resol)
                yn_1 <= yn_aux;     
                
                
                            
            END IF;     -- (reloj)
        END IF;         -- (reset - enable)
    END PROCESS proceso1;

    
    -- =======================================================================================
    --  PROCESO DE LA OTRA MITAD DE LAS OPERACIONES DE FILTRADO 
    --  (INCLUYE COEF INICIAL Y FINAL)
    -- =======================================================================================
    
    proceso2:
    PROCESS (reloj, reset)                      
    
        VARIABLE xn_a_pares : SIGNED (bits_in+1 DOWNTO 0);                          -- para sumas de pares de x[n-k] con coef simetricos
        VARIABLE REGISTRO_h : VALORES_h := (OTHERS => (OTHERS => '0'));             -- para guardar los coef (evita el registro)
        VARIABLE producto   : SIGNED (bits_in+1+bits_resol+1 DOWNTO 0);             -- prod de xn_a_pares por cada coef
        VARIABLE yn_aux     : SIGNED (bits_in+1+bits_resol+1+(num_coef-2) DOWNTO 0);-- sumas parciales
        VARIABLE signos_bi  : STD_LOGIC_VECTOR (1 DOWNTO 0);                        -- signos de los coef simetricos
        VARIABLE signos_xn  : STD_LOGIC_VECTOR (1 DOWNTO 0);                        -- Ssignos de los pares de x[n-k]
    
        -- estas variables sirven para definir cuantas operaciones el proceso para hacer las 
        -- operaciones de filtrado. El proceso comienza en el coef "indice1" y realiza operaciones 
        -- con un numero  de coeficientes dado por "margen". No se usan ctes porque solo se asignan 
        -- si num_coef supera cierto valor (sino darian error por la formula aplicada).
        VARIABLE indice1    : NATURAL;
        VARIABLE margen     : NATURAL;
        
    BEGIN
        IF (reset = '1') THEN
            yn_2 <= (OTHERS => '0');
        
        ELSIF (reloj'event AND reloj='1') THEN
            IF (enable = '1') THEN
            
                -- se usa una variable para evitar el registro de h[n], REGISTRO_h(0) = bN  
                FOR ind IN 0 TO (num_coef-1) LOOP
                    REGISTRO_h(ind) := real_a_signed (coeficiente(max_coef-num_coef+ind),bits_resol);
                END LOOP;           
    
                ----------------------------------------------------------------                                                
                -- operaciones filtrado (num_coef marca el num de prod y sumas)
                ----------------------------------------------------------------
                
                -- en este ciclo aun no se ha aplicado x[n-k] al registro de valores de x[n] 
                                
                -- como los filtros FIR son simetricos, se aprovecha esto para realizar
                -- la mitad de multiplicaciones (considerando los signos): primero se suman entre
                -- si los x[n-k] que corresponden a operaciones con coeficientes simetricos y
                -- luego se le aplica el signo (se extiende siempre el signo '0' al registro de
                -- x[n-k] porque x[n] es siempre positivo)
            
                yn_aux := (OTHERS => '0');
                
                -- operaciones con los coef simetricos (ni con el coef central ni con el inicial o final):
                
                -- (en caso de que el num de coef sea 4 o 5 las operaciones ya habrian sido realizadas
                -- por el proceso anterior)

                IF (num_coef > 6) THEN
                
                    indice1 := (num_coef/2)+(num_coef MOD 2);
                    margen  := ((num_coef-1)-(num_coef/2))/2 - (num_coef MOD 2);
                
                    FOR i IN (indice1 + margen + 1) TO (num_coef-2) LOOP
                        
                        -- signos de los bi simetricos
                        signos_bi(1):= REGISTRO_h(i)(bits_resol);
                        signos_bi(0):= REGISTRO_h(num_coef-1-i)(bits_resol);
                    
                        -- signos de los x[n-k] con coef simetricos para extension de signo
                        signos_xn(1):= REGISTRO_xn(i)(REGISTRO_xn(i)'left);
                        signos_xn(0):= REGISTRO_xn(num_coef-1-i)(REGISTRO_xn(num_coef-1-i)'left);
                    
                        CASE signos_bi IS
                            WHEN "00" =>    -- FIR tipo I - II
                                xn_a_pares:= "+"(signos_xn(1)&REGISTRO_xn(i), signos_xn(0)&REGISTRO_xn(num_coef-1-i));
                            WHEN "01" =>    -- FIR tipo III - IV
                                xn_a_pares:= "+"(signos_xn(1)&REGISTRO_xn(i), -(signos_xn(0)&REGISTRO_xn(num_coef-1-i)));
                            WHEN "10" =>    -- FIR tipo III - IV
                                xn_a_pares:= "+"(-(signos_xn(1)&REGISTRO_xn(i)), signos_xn(0)&REGISTRO_xn(num_coef-1-i));
                            WHEN OTHERS =>  -- FIR tipo I - II
                                xn_a_pares:= "+"(-(signos_xn(1)&REGISTRO_xn(i)), -(signos_xn(0)&REGISTRO_xn(num_coef-1-i)));
                        END CASE;           
                
                        producto:= "*"(xn_a_pares, ABS(REGISTRO_h(i)));
                        yn_aux := "+"(yn_aux, producto);
                    END LOOP;       
                END IF;
                    
                ----------------------------------------------------------------                                                
                -- operacion con el coef inicial y final (simetricos)
                ----------------------------------------------------------------
                
                -- El coef inicial es el que multiplica a x[n] (que aun no ha 
                -- sido registrado)             
                
                signos_bi(1):= REGISTRO_h(num_coef-1)(bits_resol);
                signos_bi(0):= REGISTRO_h(0)(bits_resol);
                
                signos_xn(0):= REGISTRO_xn(0)(REGISTRO_xn(0)'left);
                
                CASE signos_bi IS
                    WHEN "00" =>    -- FIR tipo I - II
                        xn_a_pares:= "+"("00"&SIGNED(xn), signos_xn(0)&REGISTRO_xn(0));
                    WHEN "01" =>    -- FIR tipo III - IV
                        xn_a_pares:= "+"("00"&SIGNED(xn), -(signos_xn(0)&REGISTRO_xn(0)));
                    WHEN "10" =>    -- FIR tipo III - IV
                        xn_a_pares:= "+"(-SIGNED('0'&xn), signos_xn(0)&REGISTRO_xn(0));
                    WHEN OTHERS =>  -- FIR tipo I - II
                        xn_a_pares:= "+"(-SIGNED('0'&xn), -(signos_xn(0)&REGISTRO_xn(0)));
                END CASE;                                   
                                        
                producto:= "*"(xn_a_pares, ABS(REGISTRO_h(num_coef-1)));
                yn_aux := "+"(yn_aux, producto);        

                
                --------------------------------------------------------------- 
                -- SALIDA FRACCIONAL SIN CONTROL DE SATURACION
                ---------------------------------------------------------------

                -- yn_2 agrupa los bits de la parte entera y decimal (bits_resol)
                yn_2 <= yn_aux;     
                
                
                            
            END IF;     -- (reloj)
        END IF;         -- (reset - enable)
    END PROCESS proceso2;
    
    
    -- =======================================================================================
    --  PROCESO QUE IMPLEMENTA LOS DESPLAZADORES x[n-k]
    -- =======================================================================================
        
    proceso_reg:
    PROCESS (reloj, reset)  

    BEGIN
    
        IF (reset = '1') THEN
            REGISTRO_xn <= (OTHERS => (OTHERS => '0'));
                    
        ELSIF (reloj'event AND reloj='1') THEN
            IF (enable = '1')THEN
            
                ---------------------------------------------------------------
                -- x[n-k], no se aplica al instante (agendado)  
                ---------------------------------------------------------------
                            
                -- el valor actual de x[n] entra al registro de desplazamiento por 
                -- la parte alta y el valor mas antiguo de x[n] desaparece
                
                FOR i IN 0 TO (num_coef-3) LOOP		-- (ignorado si num_coef=2)
                    REGISTRO_xn(i)<= REGISTRO_xn(i+1);
                END LOOP;
                
                REGISTRO_xn(num_coef-2)<= SIGNED('0' & xn);
                    
            END IF;
        END IF;
        
    END PROCESS proceso_reg;

    
    -- =======================================================================================
    --  CIRCUITO COMBINACIONAL QUE DETERMINA LA SALIDA
    -- =======================================================================================
    
    -- se combinan (suman) los resultados de los procesos que operan en paralelo (cada uno 
    -- realiza aprox. la mitad de operaciones de filtrado)
    
--  yn_fraccional   <= yn_1 + yn_2;     -- salida fraccional sin control de saturacion
    yn_aux          <= yn_1 + yn_2;     -- para usar en el control de saturacion de yn
    
    --------------------------------------------------------------- 
    -- SALIDA NATURAL CON CONTROL DE SATURACION
    ---------------------------------------------------------------

    -- yn_aux de tipo A(_,bits_resol): 
    --                  parte entera:  'LEFT a bits_resol 
    --                  parte decimal: (bits_resol-1) a 0           


    yn  <=  -- saturacion negativa (signo extendido)
	    	(OTHERS => '0') WHEN yn_aux(yn_aux'left) = '1' ELSE
	    	-- saturacion positiva
	    	(OTHERS => '1') WHEN yn_aux(bits_in + bits_resol)= '1' ELSE
	    	-- redondeo IEEE por arriba (si MSB de parte decimal = '1')
	    	CONV_STD_LOGIC_VECTOR(yn_aux(yn_aux'LEFT DOWNTO bits_resol),bits_in)+1 WHEN yn_aux(bits_resol-1) = '1' ELSE
	    	-- valor por defecto (sin alterar)
	    	CONV_STD_LOGIC_VECTOR(yn_aux(yn_aux'LEFT DOWNTO bits_resol),bits_in);


END arch_filtro_FIR;
--
--
--
--
--
--
--
--
--
--
--
--=============================================================================
--============ INSTANCIACION Y DECLARACION COMPONENTE filtro_FIR ==============
--=============================================================================

--COMPONENT filtro_FIR
--  GENERIC
--  (
--      bits_in     : NATURAL;      -- = bits de entrada = bits de salida
--      num_coef    : NATURAL;      -- = orden filtro + 1,  min 3 coef
--      bits_resol  : NATURAL;      -- = bits parte decimal de los coef => resolucion = 2^(-bits_resol)                 
--      coeficiente : COEFICIENTES  -- coef con matlab: "filtro_fir.m" (modificar fs, orden y fc)
--  );
--  PORT
--  (
--      reset           : IN    STD_LOGIC;                             
--      reloj           : IN    STD_LOGIC;                             
--      enable          : IN    STD_LOGIC;                             
--      xn              : IN    STD_LOGIC_VECTOR(bits_in-1 DOWNTO 0);  -- Dato de entrada síncrono
--      yn              : OUT   STD_LOGIC_VECTOR(bits_in-1 DOWNTO 0)   -- Dato de salida con control de saturacion        
--  );
--
--END COMPONENT;

--=============================================================================

--nombre_componente : filtro_FIR
--  GENERIC MAP(
--      bits_in     : NATURAL       => ,    -- = bits de entrada = bits de salida
--      num_coef    : NATURAL       => ,    -- = orden filtro + 1, min 2 coef
--      bits_resol  : NATURAL       => ,    -- = bits parte decimal de los coef => resolucion = 2^(-bits_resol)                 
--      coeficiente : COEFICIENTES  =>      -- coef con matlab: "filtro_fir.m" (modificar fs, orden y fc)
--  	)
--  PORT MAP(
--      reset   => ,        --  1bit
--      reloj   => ,        --  1bit
--      enable  => ,        --  1bit
--      xn      => ,        --  "bits_in" bit, Dato de entrada síncrono
--      yn      =>          --  "bits_in" bit, Dato de salida con control de saturacion
--  	);

--=============================================================================