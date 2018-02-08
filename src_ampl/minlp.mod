

set KEYS dimen 4;

param NONE symbolic := 'NONE';
param VAR_TYPE symbolic := 'VAR_TYPE';
param QUAD symbolic := 'QUAD';
param LIN symbolic := 'LIN';
param OBJ symbolic := 'OBJ';
param LB symbolic := 'LB';
param UB symbolic := 'UB';
param CONST symbolic := 'CONST';

param RIGHT{KEYS} symbolic;
param LEFT {KEYS} symbolic; 

set SOL dimen 4;
param SOL_REAL{SOL};
param SOL_IMAG{SOL};

set KEYS_NO_NONE_34 := {(key1, key2, key3, key4) in KEYS : key3 != NONE and key4 != NONE};

set VARIABLES := setof{(VAR_TYPE, kind, name, NONE) in KEYS}name;

set VARIABLES_ORDERED ordered by ASCII := VARIABLES;
set TRIMMER := {v1 in VARIABLES_ORDERED, v2 in VARIABLES_ORDERED: ord(v1) <= ord(v2)};

set CONSTRAINTS := 
	setof{(LB, name, NONE, NONE) in KEYS}name	
	union
	setof{(UB, name, NONE, NONE) in KEYS}name;

param CTR_LB_LEFT {name in CONSTRAINTS} := 0 + sum{(LB, name, none1, none2) in KEYS}LEFT [LB, name, none1, none2];
param CTR_LB_RIGHT{name in CONSTRAINTS} := 0 + sum{(LB, name, none1, none2) in KEYS}RIGHT[LB, name, none1, none2];
param CTR_UB_LEFT {name in CONSTRAINTS} := 0 + sum{(UB, name, none1, none2) in KEYS}LEFT [UB, name, none1, none2];
param CTR_UB_RIGHT{name in CONSTRAINTS} := 0 + sum{(UB, name, none1, none2) in KEYS}RIGHT[UB, name, none1, none2];

set CONSTRAINTS_LB := setof{(LB, name, NONE, NONE) in KEYS}name;
set CONSTRAINTS_UB := setof{(UB, name, NONE, NONE) in KEYS}name;  


set VAR_REAL_IMAG    := setof{name in VARIABLES              }(name, name&'_real', name&'_imag');
set OBJCTR_REAL_IMAG := setof{name in CONSTRAINTS union {OBJ}}(name, name&'_real', name&'_imag');

set REAL_VARIABLES   := union{name in VARIABLES,   (name, real, imag) in VAR_REAL_IMAG}   {real, imag};
set REAL_CONSTRAINTS := union{name in CONSTRAINTS, (name, real, imag) in OBJCTR_REAL_IMAG}{real};
set IMAG_CONSTRAINTS := union{name in CONSTRAINTS, (name, real, imag) in OBJCTR_REAL_IMAG}{imag};

#   (a_jk+i.b_jk)(x_j-i.y_j)(x_k+i.y_k)
# = (a_jk+i.b_jk)([x_j.x_k+y_j.y_k]+i[x_j.y_k-y_j.x_k])
# = a_jk(x_j.x_k+y_j.y_k)-b_jk(x_j.y_k-y_j.x_k) + i [b_jk(x_j.x_k+y_j.y_k)+a_jk(x_j.y_k-y_j.x_k) ]
# =
# + a_jk.x_j.x_k
# + a_jk.y_j.y_k
# - b_jk.x_j.y_k
# + b_jk.y_j.x_k
# + i [
# + b_jk.x_j.x_k
# + b_jk.y_j.y_k
# + a_jk.x_j.y_k
# - a_jk.y_j.x_k
# ]
set REAL_INDEXES :=
	# real part 
	# + a_jk.x_j.x_k
	# + a_jk.y_j.y_k
	# - b_jk.x_j.y_k
	# + b_jk.y_j.x_k
	union{
			(QUAD, objctr, var1, var2) in KEYS_NO_NONE_34, 
			(objctr, objctr_real, objctr_imag) in OBJCTR_REAL_IMAG,
			(var1, var1_real, var1_imag) in VAR_REAL_IMAG,
			(var2, var2_real, var2_imag) in VAR_REAL_IMAG
		}(
		{(objctr_real, var1_real, var2_real)} 
		union
		{(objctr_real, var1_imag, var2_imag)}
		union
		{(objctr_real, var1_real, var2_imag)}
		union
		{(objctr_real, var1_imag, var2_real)}
		)
	union
	# imag part
	# + b_jk.x_j.x_k
	# + b_jk.y_j.y_k
	# + a_jk.x_j.y_k
	# - a_jk.y_j.x_k
	union{
			(QUAD, objctr, var1, var2) in KEYS_NO_NONE_34, 
			(objctr, objctr_real, objctr_imag) in OBJCTR_REAL_IMAG,
			(var1, var1_real, var1_imag) in VAR_REAL_IMAG,
			(var2, var2_real, var2_imag) in VAR_REAL_IMAG
		}(
		{(objctr_real, var1_real, var2_real)} 
		union
		{(objctr_real, var1_imag, var2_imag)}
		union
		{(objctr_real, var1_real, var2_imag)}
		union
		{(objctr_real, var1_imag, var2_real)}
		)
;

set REAL_VARIABLES_ORDERED ordered by ASCII := REAL_VARIABLES;

set REAL_TRIMMER := {v1 in REAL_VARIABLES_ORDERED, v2 in REAL_VARIABLES_ORDERED: ord(v1) <= ord(v2)};
	
var x{REAL_VARIABLES};

# real part 
# + a_jk.x_j.x_k
# + a_jk.y_j.y_k
# - b_jk.x_j.y_k
# + b_jk.y_j.x_k
minimize CRITERION: 
+sum{(QUAD, OBJ, var1, var2) in KEYS, (var1, var1_real, var1_imag) in VAR_REAL_IMAG, (var2, var2_real, var2_imag) in VAR_REAL_IMAG}(
		+LEFT [QUAD, OBJ, var1, var2] * (if (var1_real, var2_real) in REAL_TRIMMER then x[var1_real] * x[var2_real] else x[var2_real]* x[var1_real])
		+LEFT [QUAD, OBJ, var1, var2] * (if (var1_imag, var2_imag) in REAL_TRIMMER then x[var1_imag] * x[var2_imag] else x[var2_imag]* x[var1_imag])
		-RIGHT[QUAD, OBJ, var1, var2] * (if (var1_real, var2_imag) in REAL_TRIMMER then x[var1_real] * x[var2_imag] else x[var2_imag]* x[var1_real])
		+RIGHT[QUAD, OBJ, var1, var2] * (if (var2_imag, var2_real) in REAL_TRIMMER then x[var1_imag] * x[var2_real] else x[var2_real]* x[var1_imag])  
	)
 +(if (CONST, OBJ, NONE, NONE) in KEYS then LEFT[CONST, OBJ, NONE, NONE] else 0)
	;
	
subject to constraint_real{ctr in CONSTRAINTS, (ctr, ctr_real, ctr_imag) in OBJCTR_REAL_IMAG}:
	+(if ctr in CONSTRAINTS_LB then LEFT[LB, ctr, NONE, NONE] else -Infinity)
	<=
	+sum{(QUAD, ctr, var1, var2) in KEYS, (var1, var1_real, var1_imag) in VAR_REAL_IMAG, (var2, var2_real, var2_imag) in VAR_REAL_IMAG}(
		+LEFT [QUAD, ctr, var1, var2] * (if (var1_real, var2_real) in REAL_TRIMMER then x[var1_real] * x[var2_real] else x[var2_real]* x[var1_real])
		+LEFT [QUAD, ctr, var1, var2] * (if (var1_imag, var2_imag) in REAL_TRIMMER then x[var1_imag] * x[var2_imag] else x[var2_imag]* x[var1_imag])
		-RIGHT[QUAD, ctr, var1, var2] * (if (var1_real, var2_imag) in REAL_TRIMMER then x[var1_real] * x[var2_imag] else x[var2_imag]* x[var1_real])
		+RIGHT[QUAD, ctr, var1, var2] * (if (var2_imag, var2_real) in REAL_TRIMMER then x[var1_imag] * x[var2_real] else x[var2_real]* x[var1_imag]) 
	)
 	+(if (CONST, ctr, NONE, NONE) in KEYS then LEFT[CONST, ctr, NONE, NONE] else 0)
	<=
	+(if ctr in CONSTRAINTS_UB then LEFT[UB, ctr, NONE, NONE] else +Infinity)
	;
# imag part
# + b_jk.x_j.x_k
# + b_jk.y_j.y_k
# + a_jk.x_j.y_k
# - a_jk.y_j.x_k
subject to constraint_imag{ctr in CONSTRAINTS, (ctr, ctr_real, ctr_imag) in OBJCTR_REAL_IMAG}:
	+(if ctr in CONSTRAINTS_LB then RIGHT[LB, ctr, NONE, NONE] else -Infinity)
	<=
	+sum{(QUAD, ctr, var1, var2) in KEYS, (var1, var1_real, var1_imag) in VAR_REAL_IMAG, (var2, var2_real, var2_imag) in VAR_REAL_IMAG}(
		+RIGHT[QUAD, ctr, var1, var2] * (if (var1_real, var2_real) in REAL_TRIMMER then x[var1_real] * x[var2_real] else x[var2_real] * x[var1_real])
		+RIGHT[QUAD, ctr, var1, var2] * (if (var1_imag, var2_imag) in REAL_TRIMMER then x[var1_imag] * x[var2_imag] else x[var2_imag] * x[var1_imag])
		+LEFT [QUAD, ctr, var1, var2] * (if (var1_real, var2_imag) in REAL_TRIMMER then x[var1_real] * x[var2_imag] else x[var2_imag] * x[var1_real])
		-LEFT [QUAD, ctr, var1, var2] * (if (var1_imag, var2_real) in REAL_TRIMMER then x[var1_imag] * x[var2_real] else x[var2_real] * x[var1_imag]) 
	)
 	+(if (CONST, ctr, NONE, NONE) in KEYS then RIGHT[CONST, ctr, NONE, NONE] else 0)
	<=
	+(if ctr in CONSTRAINTS_UB then RIGHT[UB, ctr, NONE, NONE] else +Infinity)
	;
	
param DUAL_REAL{ctr in CONSTRAINTS} default 0;
param DUAL_IMAG{ctr in CONSTRAINTS} default 0;


set H_INDEX dimen 2 := union{(QUAD, objctr, var1, var2) in KEYS}{(var1, var2)};

param H_REAL{(var1, var2) in H_INDEX} := 
	+(if (QUAD, OBJ, var1, var2) in KEYS then LEFT [QUAD, OBJ, var1, var2] else 0)
	+sum{(QUAD, ctr, var1, var2) in KEYS: ctr != OBJ}-(
		DUAL_REAL[ctr]*LEFT[QUAD, ctr, var1, var2]-DUAL_IMAG[ctr]*RIGHT[QUAD, ctr, var1, var2]
	)

;

param H_IMAG{(var1, var2) in H_INDEX} :=
	+(if (QUAD, OBJ, var1, var2) in KEYS then RIGHT[QUAD, OBJ, var1, var2] else 0)
	+sum{(QUAD, ctr, var1, var2) in KEYS: ctr != OBJ}-(
		DUAL_REAL[ctr]*RIGHT[QUAD, ctr, var1, var2]+DUAL_IMAG[ctr]*LEFT[QUAD, ctr, var1, var2]
	)

;