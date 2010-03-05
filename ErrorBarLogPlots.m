(* ::Package:: *)

(*: Name : ErrorBarLogPlots` *)
(*: Title : Log scale plots of data with error bars *)
(*: Summary :
	This package extends the capability provided by the standard package
	"ErrorBarPlots`" by adding functions for plots with log scales similar
	to ErrorListPlot.
  *)
(*: Context : ErrorBarLogPlots` *)
(*: Package Version : 1.0 *)
(*: Author : Frank Rice *)
(*: Copyright : 
	Original code Copyright 1988-2007 Wolfram Research, Inc.
	Modifications to code Copyright 2007, California Institute of Technology
  *)
(*: History :
      Version 1.0 by Frank Rice
  *)
(*: Keywords : plotting, 2D graphics, error bars, log scales *)
(*: Mathematica Version : 6.0 *)


BeginPackage["ErrorBarLogPlots`",{"ErrorBarPlots`"}];

ErrorListLogPlot::usage="\!\(\*RowBox[{\"ErrorListLogPlot\", \"[\", RowBox[{\"{\", "<>
	"RowBox[{RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"1\", \"TR\"]], \",\", "<>
	"SubscriptBox[StyleBox[\"dy\", \"TI\"], StyleBox[\"1\", \"TR\"]]}], \"}\"}], \",\", "<>
	"RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"2\", \"TR\"]], \",\", "<>
	"SubscriptBox[StyleBox[\"dy\", \"TI\"], StyleBox[\"2\", \"TR\"]]}], \"}\"}], \",\", "<>
	"StyleBox[\"\[Ellipsis]\", \"TI\"]}], \"}\"}], \"]\"}]\) makes a log plot of points corresponding to a list of values "<>
	"\!\(\*RowBox[{SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"1\", \"TR\"]], \",\", \" \", "<>
	"SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"2\", \"TR\"]], \",\", \" \", StyleBox[\"\[Ellipsis]\", \"TR\"]}]\), "<>
	"with corresponding error bars. The errors have magnitudes \!\(\*RowBox[{SubscriptBox[StyleBox[\"dy\", \"TI\"], "<>
	"StyleBox[\"1\", \"TR\"]], \",\", SubscriptBox[StyleBox[\"dy\", \"TI\"], StyleBox[\"2\", \"TR\"]], \",\", "<>
	"StyleBox[\"\[Ellipsis]\", \"TR\"]}]\).\n\!\(\*RowBox[{\"ErrorListLogPlot\", \"[\", RowBox[{\"{\", RowBox[{RowBox[{\"{\", "<>
	"RowBox[{RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"x\", \"TI\"], StyleBox[\"1\", \"TR\"]], \",\", "<>
	"SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"1\", \"TR\"]]}], \"}\"}], \",\", RowBox[{\"ErrorBar\", \"[\", "<>
	"SubscriptBox[StyleBox[\"err\", \"TI\"], \"1\"], \"]\"}]}], \"}\"}], \",\", RowBox[{\"{\", RowBox[{RowBox[{\"{\", "<>
	"RowBox[{SubscriptBox[StyleBox[\"x\", \"TI\"], StyleBox[\"2\", \"TR\"]], \",\", SubscriptBox[StyleBox[\"y\", \"TI\"], "<>
	"StyleBox[\"2\", \"TR\"]]}], \"}\"}], \",\", RowBox[{\"ErrorBar\", \"[\", SubscriptBox[StyleBox[\"err\", \"TI\"], "<>
	"StyleBox[\"2\", \"TR\"]], \"]\"}]}], \"}\"}], \",\", StyleBox[\"\[Ellipsis]\", \"TR\"]}], \"}\"}], \"]\"}]\) "<>
	"makes a log plot of points with specified \!\(\*StyleBox[\"x\", \"TI\"]\) and \!\(\*StyleBox[\"y\", \"TI\"]\) "<>
	"coordinates and error magnitudes.";
ErrorListLogLinearPlot::usage="\!\(\*RowBox[{\"ErrorListLogPlot\", \"[\", RowBox[{\"{\", "<>
	"RowBox[{RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"1\", \"TR\"]], \",\", "<>
	"SubscriptBox[StyleBox[\"dy\", \"TI\"], StyleBox[\"1\", \"TR\"]]}], \"}\"}], \",\", "<>
	"RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"2\", \"TR\"]], \",\", "<>
	"SubscriptBox[StyleBox[\"dy\", \"TI\"], StyleBox[\"2\", \"TR\"]]}], \"}\"}], \",\", "<>
	"StyleBox[\"\[Ellipsis]\", \"TI\"]}], \"}\"}], \"]\"}]\) makes a log-linear plot of points corresponding to a list of values "<>
	"\!\(\*RowBox[{SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"1\", \"TR\"]], \",\", \" \", "<>
	"SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"2\", \"TR\"]], \",\", \" \", StyleBox[\"\[Ellipsis]\", \"TR\"]}]\), "<>
	"with corresponding error bars. The errors have magnitudes \!\(\*RowBox[{SubscriptBox[StyleBox[\"dy\", \"TI\"], "<>
	"StyleBox[\"1\", \"TR\"]], \",\", SubscriptBox[StyleBox[\"dy\", \"TI\"], StyleBox[\"2\", \"TR\"]], \",\", "<>
	"StyleBox[\"\[Ellipsis]\", \"TR\"]}]\).\n\!\(\*RowBox[{\"ErrorListLogPlot\", \"[\", RowBox[{\"{\", RowBox[{RowBox[{\"{\", "<>
	"RowBox[{RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"x\", \"TI\"], StyleBox[\"1\", \"TR\"]], \",\", "<>
	"SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"1\", \"TR\"]]}], \"}\"}], \",\", RowBox[{\"ErrorBar\", \"[\", "<>
	"SubscriptBox[StyleBox[\"err\", \"TI\"], \"1\"], \"]\"}]}], \"}\"}], \",\", RowBox[{\"{\", RowBox[{RowBox[{\"{\", "<>
	"RowBox[{SubscriptBox[StyleBox[\"x\", \"TI\"], StyleBox[\"2\", \"TR\"]], \",\", SubscriptBox[StyleBox[\"y\", \"TI\"], "<>
	"StyleBox[\"2\", \"TR\"]]}], \"}\"}], \",\", RowBox[{\"ErrorBar\", \"[\", SubscriptBox[StyleBox[\"err\", \"TI\"], "<>
	"StyleBox[\"2\", \"TR\"]], \"]\"}]}], \"}\"}], \",\", StyleBox[\"\[Ellipsis]\", \"TR\"]}], \"}\"}], \"]\"}]\) "<>
	"makes a log-linear plot of points with specified \!\(\*StyleBox[\"x\", \"TI\"]\) and \!\(\*StyleBox[\"y\", \"TI\"]\) "<>
	"coordinates and error magnitudes.";
ErrorListLogLogPlot::usage="\!\(\*RowBox[{\"ErrorListLogPlot\", \"[\", RowBox[{\"{\", RowBox[{RowBox[{\"{\", "<>
	"RowBox[{SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"1\", \"TR\"]], \",\", SubscriptBox[StyleBox[\"dy\", \"TI\"], "<>
	"StyleBox[\"1\", \"TR\"]]}], \"}\"}], \",\", RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"y\", \"TI\"], "<>
	"StyleBox[\"2\", \"TR\"]], \",\", SubscriptBox[StyleBox[\"dy\", \"TI\"], StyleBox[\"2\", \"TR\"]]}], \"}\"}], \",\", "<>
	"StyleBox[\"\[Ellipsis]\", \"TI\"]}], \"}\"}], \"]\"}]\) makes a log-log plot of points corresponding to a list of values "<>
	"\!\(\*RowBox[{SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"1\", \"TR\"]], \",\", \" \", "<>
	"SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"2\", \"TR\"]], \",\", \" \", StyleBox[\"\[Ellipsis]\", \"TR\"]}]\), "<>
	"with corresponding error bars. The errors have magnitudes \!\(\*RowBox[{SubscriptBox[StyleBox[\"dy\", \"TI\"], "<>
	"StyleBox[\"1\", \"TR\"]], \",\", SubscriptBox[StyleBox[\"dy\", \"TI\"], StyleBox[\"2\", \"TR\"]], \",\", "<>
	"StyleBox[\"\[Ellipsis]\", \"TR\"]}]\).\n\!\(\*RowBox[{\"ErrorListLogPlot\", \"[\", RowBox[{\"{\", RowBox[{RowBox[{\"{\", "<>
	"RowBox[{RowBox[{\"{\", RowBox[{SubscriptBox[StyleBox[\"x\", \"TI\"], StyleBox[\"1\", \"TR\"]], \",\", "<>
	"SubscriptBox[StyleBox[\"y\", \"TI\"], StyleBox[\"1\", \"TR\"]]}], \"}\"}], \",\", RowBox[{\"ErrorBar\", \"[\", "<>
	"SubscriptBox[StyleBox[\"err\", \"TI\"], \"1\"], \"]\"}]}], \"}\"}], \",\", RowBox[{\"{\", RowBox[{RowBox[{\"{\", "<>
	"RowBox[{SubscriptBox[StyleBox[\"x\", \"TI\"], StyleBox[\"2\", \"TR\"]], \",\", SubscriptBox[StyleBox[\"y\", \"TI\"], "<>
	"StyleBox[\"2\", \"TR\"]]}], \"}\"}], \",\", RowBox[{\"ErrorBar\", \"[\", SubscriptBox[StyleBox[\"err\", \"TI\"], "<>
	"StyleBox[\"2\", \"TR\"]], \"]\"}]}], \"}\"}], \",\", StyleBox[\"\[Ellipsis]\", \"TR\"]}], \"}\"}], \"]\"}]\) "<>
	"makes a log-log plot of points with specified \!\(\*StyleBox[\"x\", \"TI\"]\) and \!\(\*StyleBox[\"y\", \"TI\"]\) "<>
	"coordinates and error magnitudes.";
ErrorBarScale::usage="\!\(\*RowBox[{RowBox[{\"ErrorBarScale\", \"[\"}], RowBox[{ \"]\"}]}]\) can be supplied as a value "<>
	"for the ErrorBarFunction option of the ErrorListPlot family. It scales error bar sizes to make them more visible in "<>
	"a plot.\n\!\(\*RowBox[{\"ErrorBarFunction\[RightArrow]ErrorBarScale\", \"[\", StyleBox[\"scale\", \"TI\"], \"]\"}]\) scale factor "<>
	"\!\(\*StyleBox[\"scale\", \"TI\"]\) is applied to error bars in both the x and y directions\n"<>
	"\!\(\*RowBox[{\"ErrorBarFunction\[RightArrow]ErrorBarScale\", \"[\", RowBox[{StyleBox[\"xscale\", \"TI\"], \",\", "<>
	"StyleBox[\"yscale\", \"TI\"]}], \"]\"}]\) scale factors specified independently for each of the x and the y error bars";



Begin["`Private`"];

Options[ErrorListLogLogPlot] =
Options[ErrorListLogLinearPlot] = 
Options[ErrorListLogPlot] = 
Options[elplot] = Options[ErrorListPlot];


ErrorListLogPlot[dat_, opts:OptionsPattern[]] := elplot[{False,True},dat,opts]
ErrorListLogLinearPlot[dat_, opts:OptionsPattern[]] := elplot[{True,False},dat,opts]
ErrorListLogLogPlot[dat_, opts:OptionsPattern[]] := elplot[{True,True},dat,opts]


ErrorBarScale[xscale_?NumericQ,yscale_?NumericQ] :=
	Function[{coords,errs},ErrorBarPlots`Private`ebarfun[coords,ErrorBar[xscale errs[[1]],yscale errs[[2]]]]]
ErrorBarScale[scale_?NumericQ] :=
	Function[{coords,errs},ErrorBarPlots`Private`ebarfun[coords,ErrorBar[scale errs[[1]],scale errs[[2]]]]]


(* this is the function which does all the work of log plotting *)
elplot[{xflag_,yflag_}, dat_, opts:OptionsPattern[]] :=
	Block[{xf=xflag, yf=yflag, newdat, error, newopts, ebfunc, range, p, pr, prx, pry},

		{ebfunc, range, pr} = OptionValue[{ErrorBarFunction, DataRange, PlotRange}];
	
		If[ebfunc === Automatic, ebfunc = ErrorBarPlots`Private`ebarfun];
	
		If[range === All, 
			newdat = dat /. {a_/;VectorQ[a, (NumericQ[#]||Head[#]===PlusMinus||Head[#]===ErrorBar)&] :>
				Transpose[{Range[Length[a]], a}]},
			newdat = dat /. {	
				a_/;VectorQ[a, MatchQ[#1, {_?NumericQ, _?NumericQ}] &] :>
					MapIndexed[Prepend[#,First[#2]]&, a],
				a_/;(VectorQ[a] && Length[a]>3) :> Transpose[{Range[Length[a]], a}]
				}	
			]; 
	
		newdat = newdat /. {
			{x_?NumericQ, y_?NumericQ, e_?NumericQ} :> dologs[x,y,ErrorBar[{0,0},{-e,e}]],
			{{x_,y_}, e_ErrorBar} :> dologs[x,y,makeError[e]],
			{x_/;Head[x]=!=List,y_/;Head[y]=!=List} :> handlePlusMinus[{x,y}]
			};

		newdat = newdat //. {a___,{},b___} :> {a,b};
	
		newopts = FilterRules[{opts}, Options[ListPlot]];	
			
		p = Switch[{xf,yf},
			{False,True}, ListLogPlot[newdat, newopts],
			{True,False}, ListLogLinearPlot[newdat, newopts],
			{True,True},  ListLogLogPlot[newdat, newopts]
			];

		(*Needed to handle points introduced by clipping when Joined->True*)
		error[_] := ErrorBar[{0,0},{0,0}];

		p = p /. {
			g_GraphicsComplex :> markErrors[g, ebfunc],
			l_Line :> markErrors[l,ebfunc],
			p_Point :> markErrors[p, ebfunc],
			i_Inset :> markErrors[i, ebfunc]
		};

		(*Now to fix up plot range if Automatic, All, or Full is needed *)
		If[ pr === Automatic || pr === All || pr === Full,
			{prx,pry} = PlotRange/.AbsoluteOptions[Show[p,PlotRange->pr],PlotRange];
			p = Switch[{xf,yf},
				{False,True}, ListLogPlot[newdat, PlotRange->{prx,Exp[pry]}, newopts],
				{True,False}, ListLogLinearPlot[newdat, PlotRange->{Exp[prx],pry}, newopts],
				{True,True},  ListLogLogPlot[newdat, PlotRange->{Exp[prx],Exp[pry]}, newopts]
				];
			p = p /. {
				g_GraphicsComplex :> markErrors[g, ebfunc],
				l_Line :> markErrors[l,ebfunc],
				p_Point :> markErrors[p, ebfunc],
				i_Inset :> markErrors[i, ebfunc]
			}
		];

		p
	]

makeError[ErrorBar[y_]] := ErrorBar[{0,0}, eb[y]]
makeError[ErrorBar[x_,y_]] := ErrorBar[eb[x],eb[y]]
eb[n_?Positive] := {-n,n}
eb[{n_?NumericQ, p_?NumericQ}] := {n,p}
eb[_]:={0,0}

handlePlusMinus[{x_?NumericQ, y_?NumericQ}] := 
	dologs[x,y,ErrorBar[{0,0},{0,0}]]
handlePlusMinus[{PlusMinus[x_,e_], y_?NumericQ}] := 
	dologs[x,y,ErrorBar[{-e,e},{0,0}]]
handlePlusMinus[{PlusMinus[x_,ex_], PlusMinus[y_,ey_]}] := 
	dologs[x,y,ErrorBar[{-ex,ex},{-ey,ey}]]
handlePlusMinus[{x_?NumericQ, PlusMinus[y_,ey_]}] := 
	dologs[x,y,ErrorBar[{0,0},{-ey,ey}]]
handlePlusMinus[a_] := a

dologs[x_,y_,ErrorBar[ex_,ey_]]:=
	Block[{newx=x,newy=y,newxe=ex,newye=ey,ok=True},
		If[xf, If[x>0, newx=Log[x]; newxe=ex/x, ok=False]];
		If[yf, If[y>0, newy=Log[y]; newye=ey/y, ok=False]];
		If[ok, 
			(error[N[{newx,newy}]] = ErrorBar[newxe,newye]); {x,y},
			{}
		]
	]

markErrors[GraphicsComplex[pts_, prims_, opts___], ebfunc_] := 
	GraphicsComplex[pts, prims /. {
		Line[l:{__Integer}] :> {Line[l], ebfunc[pts[[#]], error[pts[[#]]]]& /@ l},
		Line[l:{{__Integer}..}] :> {Line[l], ebfunc[pts[[#]], error[pts[[#]]]]& /@ Flatten[l]},
		Point[l:{__Integer}] :> {Point[l], ebfunc[pts[[#]], error[pts[[#]]]]& /@ l},
		Point[l:{{__Integer}..}] :> {Point[l], ebfunc[pts[[#]], error[pts[[#]]]]& /@ Flatten[l]},
		(l:Inset[obj_, pos_, a___]) :> {l, ebfunc[pts[[pos]], error[pts[[pos]]]]}		
	}, opts]

markErrors[l_Line, ebfunc_] := 
	{l, ebfunc[#, error[#]]& /@ Cases[l, {_?NumericQ, _?NumericQ}, Infinity]}

markErrors[l_Point, ebfunc_] := 
	{l, ebfunc[#, error[#]]& /@ Cases[l, {_?NumericQ, _?NumericQ}, Infinity]}

markErrors[l:Inset[obj_, pos_, a___], ebfunc_] := 
	{l, ebfunc[pos, error[pos]]}
                
End[] (* Private *) ;

EndPackage[];
