'''''''''''''''''''''''''''''''''''''''''''''''''
''''''''''''''''StatFact''''''''''''''' 	Version 1.0 - Last Revised 3rd November 2014
'''''''''''''''''''''''''''''''''''''''''''''''''

logmode +addin

'Setting up the add-in with inputs

%type = @getthistype
if %type = "NONE" then
	@uiprompt("No object found, please open a group object.")
	stop
endif
%wfname = @wfname
if %wfname="" then
	@uiprompt("No workfile is open")
	stop
endif
if @ispanel = 1 then
	@uiprompt("The procedure is not designed for panel data")
	stop
endif
%inputgroup = _this.@name
if %inputgroup = "" then  'this can occur if the group is untitled.
	%inputgroup = "_this"
endif
%numberoffactors=@getnextname("numberoffactors")	'all objects are created with the @getnextname command.
string {%numberoffactors} = "Please enter the maximum number of factors here:"
%nameofgroup=@getnextname("nameofgroup")
string {%nameofgroup} = "Please enter the name of your group here:"
%nameoffactors=@getnextname("nameoffactors")
string {%nameoffactors} = "Name of factor outputs"
!meanzerounitvariance= @hasoption("m")
!firstdifferences=@hasoption("d")
!loadingout=@hasoption("l")
!normalize=@hasoption("n")
!outputipc=@hasoption("ipc")
!outputgraphics=@hasoption("g")
if @len(%args)>0 then  'i.e. arguments were passed in, so don't put up dialog.
	%factornames = @word(%args,1)
	%fact_max=@getnextname("fact_max")
	scalar {%fact_max} = 6
	{%fact_max} = @round(@val(@word(%args,2))) 'it's rounded as we require integer values only
else
	%fact_max1=@getnextname("fact_max1")
	string {%fact_max1} = "6"
	%factornames = "stationaryfactors"
 	!result = @uidialog("Caption", "StatFact Menu", "Edit", %inputgroup, {%nameofgroup}, 21,"Edit", {%fact_max1}, {%numberoffactors},"Edit",%factornames, {%nameoffactors} , 21,"check", !firstdifferences, "Take first differences?","check", !meanzerounitvariance , "Mean Zero and Unit Variance?", "check", !normalize, "Normalize Factors (by Loadings(k_max,1))?","check", !loadingout, "Output Loadings?","check", !outputipc, "Bai and Ng (2002) Criteria?","check", !outputgraphics, "Output Graphics?")
	%fact_max=@getnextname("fact_max")
	scalar  {%fact_max}= @round(@val({%fact_max1}))
	delete {%fact_max1}
	if !result = -1 then
		delete {%nameoffactors} {%nameofgroup} {%numberoffactors} {%fact_max}
		stop
	endif
endif
if {%fact_max} >0 then ' we need a number of factors which is greater than zero
	else
	@uiprompt("Error: The maximum number of factors must be greater than 0")
	delete {%nameoffactors} {%nameofgroup} {%numberoffactors} {%fact_max}
	stop 
endif 
%temporary = %factornames+"_1"
if @isobject(%temporary) then
	@uiprompt("Error: This name for a factor object already exists in your workfile")
	delete {%nameoffactors} {%nameofgroup} {%numberoffactors} {%fact_max}
	stop 
endif
if {%fact_max} > {%inputgroup}.@count then
		@uiprompt("Error: The maximum number of factors must be less than the number of series in the group")
	delete {%nameoffactors} {%nameofgroup} {%numberoffactors} {%fact_max}
	stop 
endif
for !z = 1 to {%inputgroup}.@count
	%name = {%inputgroup}.@seriesname(!z)
	if @nas({%name})>0 then
		@uiprompt("Error: NA values not allowed in this routine")
		delete {%nameoffactors} {%nameofgroup} {%numberoffactors} {%fact_max}
		stop
	endif
next

'Begin the extraction of static factors

matrix {%inputgroup}_m = @convert({%inputgroup}) 'convert the group into a matrix
if !firstdifferences = 1 then 
	!T=@rows({%inputgroup}_m)-1
else
	!T=@rows({%inputgroup}_m)
endif
!N=@columns({%inputgroup}_m)
if !firstdifferences = 1 then 'take first differences if specified
	%inputgroup_m=@getnextname("inputgroup_m")
	matrix {%inputgroup_m} = (@subextract({%inputgroup}_m,2,1,!T+1,!N))-(@subextract({%inputgroup}_m,1,1,!T,!N))
	else
	%inputgroup_m=@getnextname("inputgroup_m")
	matrix {%inputgroup_m}={%inputgroup}_m
endif
if !meanzerounitvariance = 1 then 'standardize if specified
	for !p = 1 to @columns({%inputgroup_m})
		for !q = 1 to @rows({%inputgroup_m})
			{%inputgroup_m}(!q,!p) = ({%inputgroup_m}(!q,!p)-@mean({%inputgroup_m}.@col(!p)))/@stdev({%inputgroup_m}.@col(!p))
		next
	next
endif
delete {%inputgroup}_m 'clean up as we go along
%XX=@getnextname("XX")
if !N>=!T then 'The general structure of this part of the code follows Igor Mastens code from Banerjee et al. (2014) but is fairly standard.
	sym {%XX} = {%inputgroup_m}*@transpose({%inputgroup_m}) ' XX'
	!il = !T-{%fact_max}+1
	!iu = !T
	%va=@getnextname("va")
	matrix {%va} = @eigenvectors({%XX}) 
	%ve=@getnextname("ve")
	vector {%ve} = @eigenvalues({%XX})
	%factors=@getnextname("factors")
	matrix {%factors} = !T^(0.5)*@subextract({%va},1,!il,!T,!iu)	
	%loadings=@getnextname("loadings")
	matrix {%loadings} = (@transpose({%factors})*{%inputgroup_m})/!T
else
	sym {%XX} =(@transpose({%inputgroup_m}))*{%inputgroup_m} 'X'X
	!il = !N-{%fact_max}+1
	!iu = !N
	%va=@getnextname("va")
	matrix {%va} = @eigenvectors({%XX}) 
	%ve=@getnextname("ve")
	vector {%ve} = @eigenvalues({%XX})
	%loadings=@getnextname("loadings")
	matrix {%loadings} = (!N^(0.5))*@transpose(@subextract({%va},1,!il,!N,!iu))
	%factors=@getnextname("factors")
	matrix {%factors} = ({%inputgroup_m}*@transpose({%loadings}))/!N
endif	
if !normalize = 1 then 'Normalize as per the code of Masten
	{%factors} = {%factors}*{%loadings}({%fact_max},1)    
	{%loadings} = {%loadings}/{%loadings}({%fact_max},1)
endif
if !firstdifferences = 1 then
	smpl @first+1 @last 'account for the fact that we lose one observation
else 
	smpl @first @last
endif
%factorgroupout=@getnextname("factorgroupout")
mtos({%factors},{%factorgroupout})
for !i = 1 to {%fact_max}
	!temp = {%fact_max}-!i+1		'Renaming series for simplicity - useful later on.
	if !temp<10 then
		rename ser0!temp {%factornames}_!i
	else
		rename ser!temp {%factornames}_!i
	endif
next
if !loadingout =1 then
	if @isobject("loading_1") then
		@uiprompt("Error: You already have loading vectors in your workfile")
		delete {%va} {%XX} {%factorgroupout} {%nameoffactors} {%factors} {%ve} {%fact_max} {%nameofgroup} {%numberoffactors} {%inputgroup_m} {%loadings}
		stop 
	endif
	for !i = 1 to {%fact_max}
		vector loading_!i = @subextract({%loadings},{%fact_max}-!i+1,1,{%fact_max}-!i+1,!N)
	next
endif
delete {%va} {%XX} {%factorgroupout} 'clean up as we go along

'Begin Bai and Ng (2002) outputs 

if !outputipc = 1 then	'the format of this section is slightly different to the NgFactors.m but gives identical results
	!T=@rows({%inputgroup_m})
	!N=@columns({%inputgroup_m})
	if !N<!T then
		!Z=!N	'Corresponds to C^2 from Bai and Ng (2002
	else
		!Z=!T
	endif
	%errors = @getnextname("errors")
	matrix {%errors} = {%inputgroup_m}-({%factors}*{%loadings})
	!V_k_max=@csum(@csum(@emult({%errors},{%errors})))/(!N*!T)	'// V(k_max)
	!PC1_k_max=!V_k_max+{%fact_max}*!V_k_max*((!N+!T)/(!N*!T))*log((!N*!T)/(!N+!T))
	!PC2_k_max=!V_k_max+{%fact_max}*!V_k_max*((!N+!T)/(!N*!T))*log(!Z)
	!PC3_k_max=!V_k_max+{%fact_max}*!V_k_max*(log(!Z)/!Z)
	!IPC1_k_max=log(!V_k_max)+{%fact_max}*((!N+!T)/(!N*!T))*log((!N*!T)/(!N+!T))
	!IPC2_k_max=log(!V_k_max)+{%fact_max}*((!N+!T)/(!N*!T))*log(!Z)
	!IPC3_k_max=log(!V_k_max)+{%fact_max}*(log(!Z)/!Z)
	!AIC3_k_max=!V_k_max+{%fact_max}*!V_k_max*(2*((!N+!T-{%fact_max})/(!N*!T)))
	!BIC3_k_max=!V_k_max+{%fact_max}*!V_k_max*(((!N+!T-{%fact_max})*log(!N*!T))/(!N*!T))
	!PC1_opt = !PC1_k_max
	!PC2_opt = !PC2_k_max
	!PC3_opt = !PC3_k_max
	!IPC1_opt = !IPC1_k_max
	!IPC2_opt = !IPC2_k_max
	!IPC3_opt = !IPC3_k_max
	!AIC3_opt = !AIC3_k_max
	!BIC3_opt = !BIC3_k_max
	!k1_estim = {%fact_max}
	!k2_estim = {%fact_max}
	!k3_estim = {%fact_max}
	!Ik1_estim = {%fact_max}
	!Ik2_estim = {%fact_max}
	!Ik3_estim = {%fact_max}
	!aic3k_estim = {%fact_max}
	!bic3k_estim = {%fact_max}
	for !k = {%fact_max}-1 to 1 step -1
		{%errors} = {%inputgroup_m}-(@subextract({%factors},1,{%fact_max}-!k+1,!T,{%fact_max})*@subextract({%loadings},{%fact_max}-!k+1,1,{%fact_max},!N))
		!V_{!k} = @csum(@csum(@emult({%errors},{%errors})))/(!N*!T)
		!PC1_{!k}=!V_{!k}+!k*!V_k_max*((!N+!T)/(!N*!T))*log((!N*!T)/(!N+!T))
		!PC2_{!k}=!V_{!k}+!k*!V_k_max*((!N+!T)/(!N*!T))*log(!Z)
		!PC3_{!k}=!V_{!k}+!k*!V_k_max*(log(!Z)/!Z)
		!IPC1_{!k}=log(!V_{!k})+!k*((!N+!T)/(!N*!T))*log((!N*!T)/(!N+!T))
		!IPC2_{!k}=log(!V_{!k})+!k*((!N+!T)/(!N*!T))*log(!Z)
		!IPC3_{!k}=log(!V_{!k})+!k*(log(!Z)/!Z)
		!AIC3_{!k}=!V_{!k}+!k*!V_k_max*(2*((!N+!T-!k)/(!N*!T)))
		!BIC3_{!k}=!V_{!k}+!k*!V_k_max*(((!N+!T-!k)*log(!N*!T))/(!N*!T))
		if !PC1_{!k}<!PC1_opt then
			!k1_estim= !k
			!PC1_opt=!PC1_{!k}
		endif
		if !PC2_{!k}<!PC2_opt then
			!k2_estim= !k
			!PC2_opt=!PC2_{!k}
		endif
		if !PC3_{!k}<!PC3_opt then
			!k3_estim= !k
			!PC3_opt=!PC3_{!k}
		endif
		if !IPC1_{!k}<!IPC1_opt then
			!Ik1_estim= !k
			!IPC1_opt=!IPC1_{!k}
		endif
		if !IPC2_{!k}<!IPC2_opt then
			!Ik2_estim= !k
			!IPC2_opt=!IPC2_{!k}
		endif
		if !IPC3_{!k}<!IPC3_opt then
			!Ik3_estim= !k
			!IPC3_opt=!IPC3_{!k}
		endif
		if !AIC3_{!k}<!AIC3_opt then
			!aic3k_estim= !k
			!AIC3_opt=!AIC3_{!k}
		endif
		if !BIC3_{!k}<!BIC3_opt then
			!bic3k_estim= !k
			!BIC3_opt=!BIC3_{!k}
		endif
	next
	%baingresults=@getnextname("baingresults")
	table {%baingresults} 'create and fill the output table for this section. Other output structures entirely possible.
	{%baingresults}(1,1) = "Bai and Ng (2002) Results"
	{%baingresults}(2,1) = "-----------------------------------"
	{%baingresults}(4,2) = "PC1"
	{%baingresults}(4,3) = "PC2"
	{%baingresults}(4,4) = "PC3"
	{%baingresults}(4,5) = "IPC1"
	{%baingresults}(4,6) = "IPC2"
	{%baingresults}(4,7) = "IPC3"
	{%baingresults}(4,8) = "AIC3"
	{%baingresults}(4,9) = "BIC3"
	{%baingresults}(5,1) = {%fact_max}
	{%baingresults}(5,2) = !PC1_k_max
	{%baingresults}(5,3) = !PC2_k_max
	{%baingresults}(5,4) = !PC3_k_max
	{%baingresults}(5,5) = !IPC1_k_max
	{%baingresults}(5,6) = !IPC2_k_max
	{%baingresults}(5,7) = !IPC3_k_max
	{%baingresults}(5,8) = !AIC3_k_max
	{%baingresults}(5,9) = !BIC3_k_max
	for !k = {%fact_max}-1 to 1 step -1
		{%baingresults}(5+{%fact_max}-!k,1) = !k
		{%baingresults}(5+{%fact_max}-!k,2)= !PC1_{!k}
		{%baingresults}(5+{%fact_max}-!k,3)= !PC2_{!k}	
		{%baingresults}(5+{%fact_max}-!k,4)= !PC3_{!k}
		{%baingresults}(5+{%fact_max}-!k,5)= !IPC1_{!k}
		{%baingresults}(5+{%fact_max}-!k,6)= !IPC2_{!k}
		{%baingresults}(5+{%fact_max}-!k,7)= !IPC3_{!k}
		{%baingresults}(5+{%fact_max}-!k,8)= !AIC3_{!k}
		{%baingresults}(5+{%fact_max}-!k,9)= !BIC3_{!k}	
	next
	{%baingresults}(5+{%fact_max},2) = "--------"
	{%baingresults}(5+{%fact_max},3) = "--------"
	{%baingresults}(5+{%fact_max},4) = "--------"
	{%baingresults}(5+{%fact_max},5) = "--------"
	{%baingresults}(5+{%fact_max},6) = "--------"
	{%baingresults}(5+{%fact_max},7) = "--------"
	{%baingresults}(5+{%fact_max},8) = "--------"
	{%baingresults}(5+{%fact_max},9) = "--------"
	{%baingresults}(5+{%fact_max}+1,1) = "Optimal"
	{%baingresults}(5+{%fact_max}+1,2) = !k1_estim
	{%baingresults}(5+{%fact_max}+1,3) = !k2_estim
	{%baingresults}(5+{%fact_max}+1,4) = !k3_estim
	{%baingresults}(5+{%fact_max}+1,5) = !Ik1_estim
	{%baingresults}(5+{%fact_max}+1,6) = !Ik2_estim
	{%baingresults}(5+{%fact_max}+1,7) = !Ik3_estim
	{%baingresults}(5+{%fact_max}+1,8) = !aic3k_estim
	{%baingresults}(5+{%fact_max}+1,9) = !bic3k_estim
	delete {%errors}
	show {%baingresults}
endif

if !outputgraphics = 1 then
	%factor_group = @getnextname("factor_group")
	group {%factor_group}
	for !k = 1 to {%fact_max}
		{%factor_group}.add {%factornames}_!k
	next
	%factor_graph = @getnextname("factor_graph")
	graph {%factor_graph}.line(m) {%factor_group}
	{%factor_graph}.options fillcolor(white) backcolor(white) framecolor(black) fillfade(none) backfade(none) 
	{%factor_graph}.datelabel format(yyyy:mm)
	{%factor_graph}.setelem(1) linecolor(blue)
	{%factor_graph}.addtext(t, font(arial, 20)) "Factors"
	matrix {%loadings}trans=@transpose({%loadings})
	%loadings_graph =@getnextname("loadings_graph")
	freeze({%loadings_graph}) {%loadings}trans.dot
	for !k = 1 to {%fact_max}
		{%loadings_graph}.setelem(!k) legend(Loading to Factor !k)
	next
	{%loadings_graph}.addtext(t, font(arial, 20)) "Loadings"
	{%loadings_graph}.options fillcolor(white) backcolor(white) framecolor(black) fillfade(none) backfade(none) 
	%stationary_screeplot = @getnextname("stationary_screeplot")
	vector(@rows({%ve})) {%stationary_screeplot}
	for !t = 1 to @rows({%ve})
		{%stationary_screeplot}(!t) = {%ve}(@rows({%ve})-!t+1)
	next
	%stationary_screeplotg = @getnextname("stationary_screeplotg")
	freeze({%stationary_screeplotg}) {%stationary_screeplot}.line
	{%stationary_screeplotg}.addtext(t, font(arial, 20)) "Scree Plot: Stationary Variable"
	{%stationary_screeplotg}.setelem(1) linecolor(blue) legend(Eigenvalue Number)
	{%stationary_screeplotg}.legend font("arial",18)
	{%stationary_screeplotg}.setelem(1) legend(Eigenvalue Number)
	{%stationary_screeplotg}.options fillcolor(white) backcolor(white) framecolor(black) fillfade(none) backfade(none) 
	{%stationary_screeplotg}.align(2,1.5,1.5)
	{%stationary_screeplotg}.legend +display
	!sum_eig = 0		
	%pc_expl_stat=@getnextname("pc_expl_stat")										
	vector({%fact_max}) {%pc_expl_stat}
	for !r = 1 to @rows({%ve})
		!sum_eig = {%ve}(@rows({%ve})-!r+1)+!sum_eig
	next
	for !k = 1 to {%fact_max}
		{%pc_expl_stat}(!k)=({%ve}(@rows({%ve})-!k+1)/!sum_eig)*100
	next
	%pc_explained_stat_g=@getnextname("pc_explained_stat_g")
	freeze({%pc_explained_stat_g}) {%pc_expl_stat}.bar
	{%pc_explained_stat_g}.setelem(1) linecolor(blue)
	{%pc_explained_stat_g}.legend font("arial",18)
	{%pc_explained_stat_g}.axis(l) format(suffix=%)
	{%pc_explained_stat_g}.addtext(t, font(arial, 20)) "Percent Explained: Stationary Factors"
	{%pc_explained_stat_g}.setelem(1) legend(Factor Number)
	{%pc_explained_stat_g}.legend +display
	{%pc_explained_stat_g}.options fillcolor(white) backcolor(white) framecolor(black) fillfade(none) backfade(none) 
	{%pc_explained_stat_g}.align(2,1.5,1.5)
	show {%pc_explained_stat_g}
	show {%stationary_screeplotg}
	show {%factor_graph}
	show {%loadings_graph}
	delete {%stationary_screeplot} {%pc_expl_stat} {%factor_group} {%loadings}trans
endif 
smpl @all
delete {%nameoffactors} {%factors} {%ve} {%fact_max} {%nameofgroup} {%numberoffactors} {%inputgroup_m} {%loadings}


