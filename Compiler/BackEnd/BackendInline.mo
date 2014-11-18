/*
 * This file is part of OpenModelica.
 *
 * Copyright (c) 1998-2014, Open Source Modelica Consortium (OSMC),
 * c/o Linköpings universitet, Department of Computer and Information Science,
 * SE-58183 Linköping, Sweden.
 *
 * All rights reserved.
 *
 * THIS PROGRAM IS PROVIDED UNDER THE TERMS OF GPL VERSION 3 LICENSE OR
 * THIS OSMC PUBLIC LICENSE (OSMC-PL) VERSION 1.2.
 * ANY USE, REPRODUCTION OR DISTRIBUTION OF THIS PROGRAM CONSTITUTES
 * RECIPIENT'S ACCEPTANCE OF THE OSMC PUBLIC LICENSE OR THE GPL VERSION 3,
 * ACCORDING TO RECIPIENTS CHOICE.
 *
 * The OpenModelica software and the Open Source Modelica
 * Consortium (OSMC) Public License (OSMC-PL) are obtained
 * from OSMC, either from the above address,
 * from the URLs: http://www.ida.liu.se/projects/OpenModelica or
 * http://www.openmodelica.org, and in the OpenModelica distribution.
 * GNU version 3 is obtained from: http://www.gnu.org/copyleft/gpl.html.
 *
 * This program is distributed WITHOUT ANY WARRANTY; without
 * even the implied warranty of  MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE, EXCEPT AS EXPRESSLY SET FORTH
 * IN THE BY RECIPIENT SELECTED SUBSIDIARY LICENSE CONDITIONS OF OSMC-PL.
 *
 * See the full OSMC Public License conditions for more details.
 *
 */

encapsulated package BackendInline
" file:        BackendInline.mo
  package:     BackendInline
  description: inline functions

  RCS: $Id$

  This module contains data structures and functions for inline functions.

  The entry point is the inlineCalls function, or inlineCallsInFunctions
  "

public import Absyn;
public import BackendDAE;
public import BaseHashTable;
public import DAE;
public import FCore;
public import HashTableCG;
public import Inline;
public import SCode;
public import Values;

protected import Ceval;
protected import ClassInf;
protected import ComponentReference;
protected import Config;
protected import Debug;
protected import Error;
protected import Expression;
protected import ExpressionDump;
protected import ExpressionSimplify;
protected import Flags;
protected import List;
protected import Types;
protected import VarTransform;

// =============================================================================
// late inline functions stuff
//
// =============================================================================
public function lateInlineFunction
  input BackendDAE.BackendDAE inDAE;
  output BackendDAE.BackendDAE outDAE;
algorithm
  outDAE := inlineCalls({DAE.NORM_INLINE(), DAE.AFTER_INDEX_RED_INLINE()}, inDAE);
end lateInlineFunction;

// =============================================================================
// inline calls stuff
//
// =============================================================================

public function inlineCalls
"searches for calls where the inline flag is true, and inlines them"
  input list<DAE.InlineType> inITLst;
  input BackendDAE.BackendDAE inBackendDAE;
  output BackendDAE.BackendDAE outBackendDAE;
algorithm
  outBackendDAE := matchcontinue(inBackendDAE)
    local
      list<DAE.InlineType> itlst;
      Inline.Functiontuple tpl;
      BackendDAE.EqSystems eqs;
      BackendDAE.Shared shared;

    case BackendDAE.DAE(eqs, shared as BackendDAE.SHARED())
      algorithm
        tpl := (SOME(shared.functionTree), inITLst);
        eqs := List.map1(eqs, inlineEquationSystem, tpl);
        shared.knownVars := inlineVariables(shared.knownVars, tpl);
        shared.externalObjects := inlineVariables(shared.externalObjects, tpl);
        shared.initialEqs := inlineEquationArray(shared.initialEqs, tpl);
        shared.removedEqs := inlineEquationArray(shared.removedEqs, tpl);
        shared.eventInfo := inlineEventInfo(shared.eventInfo, tpl);
      then
        BackendDAE.DAE(eqs, shared);

    else
      algorithm
        true := Flags.isSet(Flags.FAILTRACE);
        Debug.traceln("Inline.inlineCalls failed");
      then
        fail();
  end matchcontinue;
end inlineCalls;

protected function inlineEquationSystem
  input BackendDAE.EqSystem eqs;
  input Inline.Functiontuple tpl;
  output BackendDAE.EqSystem oeqs;
algorithm
  oeqs := match (eqs,tpl)
    local
      BackendDAE.EqSystem syst;
      BackendDAE.Variables orderedVars;
      BackendDAE.EquationArray orderedEqs;
      BackendDAE.Matching matching;
      Boolean b1,b2;
      BackendDAE.StateSets stateSets;
      BackendDAE.BaseClockPartitionKind partitionKind;

    case (syst as BackendDAE.EQSYSTEM(orderedVars=orderedVars,orderedEqs=orderedEqs,matching=matching,stateSets=stateSets,partitionKind=partitionKind),_)
      equation
        (orderedVars,b1) = inlineVariables(orderedVars,tpl);
        (orderedEqs,b2) = inlineEquationArray(orderedEqs,tpl);
        syst = if b1 or b2 then BackendDAE.EQSYSTEM(orderedVars,orderedEqs,NONE(),NONE(),matching,stateSets,partitionKind) else syst;
      then
        syst;
  end match;
end inlineEquationSystem;

protected function inlineEquationArray "
function: inlineEquationArray
  inlines function calls in an equation array"
  input BackendDAE.EquationArray inEquationArray;
  input Inline.Functiontuple inElementList;
  output BackendDAE.EquationArray outEquationArray;
  output Boolean oInlined;
algorithm
  (outEquationArray,oInlined) := matchcontinue(inEquationArray,inElementList)
    local
      Inline.Functiontuple fns;
      Integer i1,i2,size;
      array<Option<BackendDAE.Equation>> eqarr;
    case(BackendDAE.EQUATION_ARRAY(size,i1,i2,eqarr),fns)
      equation
        oInlined = inlineEquationOptArray(eqarr,i2,fns);
      then
        (BackendDAE.EQUATION_ARRAY(size,i1,i2,eqarr),oInlined);
    else
      equation
        true = Flags.isSet(Flags.FAILTRACE);
        Debug.trace("Inline.inlineEquationArray failed\n");
      then
        fail();
  end matchcontinue;
end inlineEquationArray;

protected function inlineEquationOptArray
"functio: inlineEquationrOptArray
  inlines calls in a equation option"
  input array<Option<BackendDAE.Equation>> inEqnArray;
  input Integer arraysize;
  input Inline.Functiontuple fns;
  output Boolean oInlined := false;
protected
  Option<BackendDAE.Equation> eqn;
  Boolean inlined;
algorithm
  for i in 1:arraysize loop
    (eqn, inlined) := inlineEqOpt(inEqnArray[i], fns);

    if inlined then
      arrayUpdate(inEqnArray, i, eqn);
      oInlined := true;
    end if;
  end for;
end inlineEquationOptArray;

protected function inlineEqOpt "
function: inlineEqOpt
  inlines function calls in equations"
  input Option<BackendDAE.Equation> inEquationOption;
  input Inline.Functiontuple inElementList;
  output Option<BackendDAE.Equation> outEquationOption;
  output Boolean inlined;
algorithm
  (outEquationOption,inlined) := match(inEquationOption,inElementList)
    local
      BackendDAE.Equation eqn;
      Boolean b;
    case(NONE(),_) then (NONE(),false);
    case(SOME(eqn),_)
      equation
        (eqn,b) = inlineEq(eqn,inElementList);
      then
        (SOME(eqn),b);
  end match;
end inlineEqOpt;

protected function inlineEq "
  inlines function calls in equations"
  input BackendDAE.Equation inEquation;
  input Inline.Functiontuple fns;
  output BackendDAE.Equation outEquation;
  output Boolean inlined;
algorithm
  (outEquation,inlined) := matchcontinue(inEquation,fns)
    local
      DAE.Exp e,e_1,e1,e1_1,e2,e2_1;
      Integer size;
      list<DAE.Exp> explst;
      DAE.ComponentRef cref;
      BackendDAE.WhenEquation weq,weq_1;
      DAE.ElementSource source;
      list<Integer> dimSize;
      DAE.Algorithm alg;
      list<DAE.Statement> stmts,stmts1,assrtLst;
      list<BackendDAE.Equation> eqns;
      list<list<BackendDAE.Equation>> eqnslst;
      Boolean b1,b2,b3;
      DAE.Expand crefExpand;
      BackendDAE.EquationAttributes attr;

    case(BackendDAE.EQUATION(e1,e2,source,attr),_)
      equation
        (e1_1,source,b1,_) = Inline.inlineExp(e1,fns,source);
        (e2_1,source,b2,_) = Inline.inlineExp(e2,fns,source);
        true = b1 or b2;
      then
       (BackendDAE.EQUATION(e1_1,e2_1,source,attr),true);

    case(BackendDAE.ARRAY_EQUATION(dimSize,e1,e2,source,attr),_)
      equation
        (e1_1,source,b1,_) = Inline.inlineExp(e1,fns,source);
        (e2_1,source,b2,_) = Inline.inlineExp(e2,fns,source);
        true = b1 or b2;
      then
        (BackendDAE.ARRAY_EQUATION(dimSize,e1_1,e2_1,source,attr),true);

    case(BackendDAE.SOLVED_EQUATION(cref,e,source,attr),_)
      equation
        (e_1,source,true,_) = Inline.inlineExp(e,fns,source);
      then
        (BackendDAE.SOLVED_EQUATION(cref,e_1,source,attr),true);

    case(BackendDAE.RESIDUAL_EQUATION(e,source,attr),_)
      equation
        (e_1,source,true,_) = Inline.inlineExp(e,fns,source);
      then
        (BackendDAE.RESIDUAL_EQUATION(e_1,source,attr),true);

    case(BackendDAE.ALGORITHM(size,alg as DAE.ALGORITHM_STMTS(statementLst=stmts),source,crefExpand,attr),_)
      equation
        (stmts1,true) = Inline.inlineStatements(stmts,fns,{},false);
        alg = DAE.ALGORITHM_STMTS(stmts1);
      then
        (BackendDAE.ALGORITHM(size,alg,source,crefExpand,attr),true);

    case(BackendDAE.WHEN_EQUATION(size,weq,source,attr),_)
      equation
        (weq_1,source,true) = inlineWhenEq(weq,fns,source);
      then
        (BackendDAE.WHEN_EQUATION(size,weq_1,source,attr),true);

    case(BackendDAE.COMPLEX_EQUATION(size,e1,e2,source,attr),_)
      equation
        (e1_1,source,b1,_) = Inline.inlineExp(e1,fns,source);
        (e2_1,source,b2,_) = Inline.inlineExp(e2,fns,source);
        true = b1 or b2;
      then
        (BackendDAE.COMPLEX_EQUATION(size,e1_1,e2_1,source,attr),true);

    case(BackendDAE.IF_EQUATION(explst,eqnslst,eqns,source,attr),_)
      equation
        (explst,source,b1) = Inline.inlineExps(explst,fns,source);
        (eqnslst,b2) = inlineEqsLst(eqnslst,fns,{},false);
        (eqns,b3) = inlineEqs(eqns,fns,{},false);
        true = b1 or b2 or b3;
      then
        (BackendDAE.IF_EQUATION(explst,eqnslst,eqns,source,attr),true);
    else
      then
        (inEquation,false);
  end matchcontinue;
end inlineEq;

protected function inlineEqsLst
  input list<list<BackendDAE.Equation>> inEqnsList;
  input Inline.Functiontuple inFunctions;
  input list<list<BackendDAE.Equation>> iAcc;
  input Boolean iInlined;
  output list<list<BackendDAE.Equation>> outEqnsList;
  output Boolean OInlined;
algorithm
  (outEqnsList,OInlined) := match(inEqnsList,inFunctions,iAcc,iInlined)
    local
      list<BackendDAE.Equation> eqn;
      list<list<BackendDAE.Equation>> rest,acc;
      Boolean inlined;
    case ({},_,_,_) then (listReverse(iAcc),iInlined);
    case (eqn::rest,_,_,_)
      equation
        (eqn,inlined) = inlineEqs(eqn,inFunctions,{},false);
        (acc,inlined) = inlineEqsLst(rest,inFunctions,eqn::iAcc,inlined or iInlined);
      then
        (acc,inlined);
  end match;
end inlineEqsLst;

public function inlineEqs
  input list<BackendDAE.Equation> inEqnsList;
  input Inline.Functiontuple inFunctions;
  input list<BackendDAE.Equation> iAcc;
  input Boolean iInlined;
  output list<BackendDAE.Equation> outEqnsList;
  output Boolean OInlined;
algorithm
  (outEqnsList,OInlined) := match(inEqnsList,inFunctions,iAcc,iInlined)
    local
      BackendDAE.Equation eqn;
      list<BackendDAE.Equation> rest,acc;
      Boolean inlined;
    case ({},_,_,_) then (listReverse(iAcc),iInlined);
    case (eqn::rest,_,_,_)
      equation
        (eqn,inlined) = inlineEq(eqn,inFunctions);
        (acc,inlined) = inlineEqs(rest,inFunctions,eqn::iAcc,inlined or iInlined);
      then
        (acc,inlined);
  end match;
end inlineEqs;

protected function inlineWhenEq
"inlines function calls in when equations"
  input BackendDAE.WhenEquation inWhenEquation;
  input Inline.Functiontuple fns;
  input DAE.ElementSource inSource;
  output BackendDAE.WhenEquation outWhenEquation;
  output DAE.ElementSource outSource;
  output Boolean inlined;
algorithm
  (outWhenEquation,outSource,inlined) := matchcontinue(inWhenEquation,fns,inSource)
    local
      DAE.ComponentRef cref;
      DAE.Exp e,e_1,cond;
      BackendDAE.WhenEquation weq,weq_1;
      DAE.ElementSource source;
      Boolean b1,b2,b3;
      list<DAE.Statement> assrtLst;
    case (BackendDAE.WHEN_EQ(cond,cref,e,NONE()),_,_)
      equation
        (e_1,source,b1,_) = Inline.inlineExp(e,fns,inSource);
        (cond,source,b2,_) = Inline.inlineExp(cond,fns,source);
        true = b1 or b2;
      then
        (BackendDAE.WHEN_EQ(cond,cref,e_1,NONE()),source,true);
    case (BackendDAE.WHEN_EQ(cond,cref,e,SOME(weq)),_,_)
      equation
        (e_1,source,b1,_) = Inline.inlineExp(e,fns,inSource);
        (cond,source,b2,_) = Inline.inlineExp(cond,fns,source);
        (weq_1,source,b3) = inlineWhenEq(weq,fns,source);
        true = b1 or b2 or b3;
      then
        (BackendDAE.WHEN_EQ(cond,cref,e_1,SOME(weq_1)),source,true);
    else
      then
        (inWhenEquation,inSource,false);
  end matchcontinue;
end inlineWhenEq;

protected function inlineVariables
"inlines function calls in variables"
  input BackendDAE.Variables inVariables;
  input Inline.Functiontuple inElementList;
  output BackendDAE.Variables outVariables;
  output Boolean inlined;
algorithm
  (outVariables,inlined) := matchcontinue(inVariables,inElementList)
    local
      Inline.Functiontuple fns;
      array<list<BackendDAE.CrefIndex>> crefind;
      Integer i1,i2,i3,i4;
      array<Option<BackendDAE.Var>> vararr;
    case(BackendDAE.VARIABLES(crefind,BackendDAE.VARIABLE_ARRAY(i3,i4,vararr),i1,i2),fns)
      equation
        inlined = inlineVarOptArray(1,vararr,i4,fns,false);
      then
        (BackendDAE.VARIABLES(crefind,BackendDAE.VARIABLE_ARRAY(i3,i4,vararr),i1,i2),inlined);
    else
      equation
        true = Flags.isSet(Flags.FAILTRACE);
        Debug.trace("Inline.inlineVariables failed\n");
      then
        fail();
  end matchcontinue;
end inlineVariables;

protected function updateArrayCond
  input Boolean cond;
  input array<Type_a> inArr;
  input Integer index;
  input Type_a value;
  replaceable type Type_a subtypeof Any;
algorithm
  _ := match(cond,inArr,index,value)
    case(false,_,_,_) then ();
    case(true,_,_,_)
      equation
        arrayUpdate(inArr,index,value);
      then
        ();
  end match;
end updateArrayCond;

protected function inlineVarOptArray
"functio: inlineVarOptArray
  inlines calls in a variable option"
  input Integer index;
  input array<Option<BackendDAE.Var>> inVarArray;
  input Integer arraysize;
  input Inline.Functiontuple fns;
  input Boolean iInlined;
  output Boolean oInlined;
algorithm
  oInlined := inlineVarOptArrayWork(index > arraysize,index,inVarArray,arraysize,fns,iInlined);
end inlineVarOptArray;

protected function inlineVarOptArrayWork
"functio: inlineVarOptArray
  inlines calls in a variable option"
  input Boolean stop;
  input Integer index;
  input array<Option<BackendDAE.Var>> inVarArray;
  input Integer arraysize;
  input Inline.Functiontuple fns;
  input Boolean iInlined;
  output Boolean oInlined;
algorithm
  oInlined := match (stop,index,inVarArray,arraysize,fns,iInlined)
    local
      Option<BackendDAE.Var> var;
      Boolean b;
    case (true,_,_,_,_,_)
      then iInlined;
    else
      equation
        var = inVarArray[index];
        (var,b) = inlineVarOpt(var,fns);
        updateArrayCond(b,inVarArray,index,var);
      then inlineVarOptArrayWork(index+1 > arraysize,index+1,inVarArray,arraysize,fns,b or iInlined);
  end match;
end inlineVarOptArrayWork;

protected function inlineVarOpt
"functio: inlineVarOpt
  inlines calls in a variable option"
  input Option<BackendDAE.Var> inVarOption;
  input Inline.Functiontuple fns;
  output Option<BackendDAE.Var> outVarOption;
  output Boolean inlined;
algorithm
  (outVarOption,inlined) := match(inVarOption,fns)
    local
      BackendDAE.Var var;
      Boolean b;
    case(NONE(),_) then (NONE(),false);
    case(SOME(var),_)
      equation
        (var,b) = inlineVar(var,fns);
      then
        (SOME(var),b);
  end match;
end inlineVarOpt;

protected function inlineVar
"functio: inlineVar
  inlines calls in a variable"
  input BackendDAE.Var inVar;
  input Inline.Functiontuple inElementList;
  output BackendDAE.Var outVar;
  output Boolean inlined;
algorithm
  (outVar,inlined) := match(inVar,inElementList)
    local
      Inline.Functiontuple fns;
      DAE.ComponentRef varName;
      BackendDAE.VarKind varKind;
      DAE.VarDirection varDirection;
      DAE.VarParallelism varParallelism;
      BackendDAE.Type varType;
      Option<Values.Value> bindValue;
      DAE.InstDims arrayDim;
      Option<DAE.VariableAttributes> values,values1;
	  Option<BackendDAE.TearingSelect> ts;
      Option<SCode.Comment> comment;
      DAE.ConnectorType ct;
      BackendDAE.Var var;
      DAE.ElementSource source;
      Option<DAE.Exp> bind;
      Boolean b1,b2;
    case(BackendDAE.VAR(varName,varKind,varDirection,varParallelism,varType,bind,bindValue,arrayDim,source,values,ts,comment,ct),fns)
      equation
        (bind,source,b1) = Inline.inlineExpOpt(bind,fns,source);
        (values1,source,b2) = Inline.inlineStartAttribute(values,source,fns);
      then
        (BackendDAE.VAR(varName,varKind,varDirection,varParallelism,varType,bind,bindValue,arrayDim,source,values1,ts,comment,ct),b1 or b2);
    case(var,_) then (var,false);
  end match;
end inlineVar;

protected function inlineEventInfo "inlines function calls in event info"
  input BackendDAE.EventInfo inEventInfo;
  input Inline.Functiontuple inElementList;
  output BackendDAE.EventInfo outEventInfo;
algorithm
  outEventInfo := matchcontinue(inEventInfo, inElementList)
    local
      Inline.Functiontuple fns;
      list<BackendDAE.WhenClause> wclst, wclst_1;
      list<BackendDAE.ZeroCrossing> zclst, zclst_1, relations, samples;
      Integer numberOfMathEvents;
      BackendDAE.EventInfo ev;
      Boolean b1, b2, b3;
      list<BackendDAE.TimeEvent> timeEvents;

    case(BackendDAE.EVENT_INFO(timeEvents, wclst, zclst, samples, relations, numberOfMathEvents), fns) equation
      (wclst_1, b1) = inlineWhenClauses(wclst, fns, {}, false);
      (zclst_1, b2) = inlineZeroCrossings(zclst, fns, {}, false);
      (relations, b3) = inlineZeroCrossings(relations, fns, {}, false);
      ev = if b1 or b2 or b3 then BackendDAE.EVENT_INFO(timeEvents, wclst_1, zclst_1, samples, relations, numberOfMathEvents) else inEventInfo;
    then ev;

    else
      equation
        true = Flags.isSet(Flags.FAILTRACE);
        Debug.trace("Inline.inlineEventInfo failed\n");
      then fail();
  end matchcontinue;
end inlineEventInfo;

protected function inlineZeroCrossings "inlines function calls in zero crossings"
  input list<BackendDAE.ZeroCrossing> inStmts;
  input Inline.Functiontuple fns;
  input list<BackendDAE.ZeroCrossing> iAcc;
  input Boolean iInlined;
  output list<BackendDAE.ZeroCrossing> outStmts;
  output Boolean oInlined;
algorithm
  (outStmts, oInlined) := match (inStmts, fns, iAcc, iInlined)
    local
      BackendDAE.ZeroCrossing zc;
      list<BackendDAE.ZeroCrossing> rest, stmts;
      Boolean b;

    case ({}, _, _, _)
    then (listReverse(iAcc), iInlined);

    case (zc::rest, _, _, _) equation
      (zc, b) = inlineZeroCrossing(zc, fns);
      (stmts, b) = inlineZeroCrossings(rest, fns, zc::iAcc, b or iInlined);
    then (stmts, b);
  end match;
end inlineZeroCrossings;

protected function inlineZeroCrossing "inlines function calls in a zero crossing"
  input BackendDAE.ZeroCrossing inZeroCrossing;
  input Inline.Functiontuple inElementList;
  output BackendDAE.ZeroCrossing outZeroCrossing;
  output Boolean oInlined;
algorithm
  (outZeroCrossing, oInlined) := matchcontinue(inZeroCrossing, inElementList)
    local
      Inline.Functiontuple fns;
      DAE.Exp e, e_1;
      list<Integer> ilst1, ilst2;
      list<DAE.Statement> assrtLst;

    case(BackendDAE.ZERO_CROSSING(e, ilst1, ilst2), fns) equation
      (e_1, _, true, _) = Inline.inlineExp(e, fns, DAE.emptyElementSource/*TODO: Propagate operation info*/);
    then (BackendDAE.ZERO_CROSSING(e_1, ilst1, ilst2), true);

    else (inZeroCrossing, false);
  end matchcontinue;
end inlineZeroCrossing;

protected function inlineWhenClauses
"inlines function calls in reinit statements"
  input list<BackendDAE.WhenClause> inStmts;
  input Inline.Functiontuple fns;
  input list<BackendDAE.WhenClause> iAcc;
  input Boolean iInlined;
  output list<BackendDAE.WhenClause> outStmts;
  output Boolean oInlined;
algorithm
  (outStmts,oInlined) := match (inStmts,fns,iAcc,iInlined)
    local
      BackendDAE.WhenClause wc;
      list<BackendDAE.WhenClause> rest,stmts;
      Boolean b;

    case ({},_,_,_) then (listReverse(iAcc),iInlined);
    case (wc::rest,_,_,_)
      equation
        (wc,b) = inlineWhenClause(wc,fns);
        (stmts,b) = inlineWhenClauses(rest,fns,wc::iAcc,b or iInlined);
      then
        (stmts,b);
  end match;
end inlineWhenClauses;

protected function inlineWhenClause
"inlines function calls in a when clause"
  input BackendDAE.WhenClause inWhenClause;
  input Inline.Functiontuple inElementList;
  output BackendDAE.WhenClause outWhenClause;
  output Boolean inlined;
algorithm
  (outWhenClause,inlined) := matchcontinue(inWhenClause,inElementList)
    local
      Inline.Functiontuple fns;
      DAE.Exp e,e_1;
      list<BackendDAE.WhenOperator> rslst,rslst_1;
      Option<Integer> io;
      Boolean b1,b2;
      list<DAE.Statement> assrtLst;
    case(BackendDAE.WHEN_CLAUSE(e,rslst,io),fns)
      equation
        (e_1,_,b1,_) = Inline.inlineExp(e,fns,DAE.emptyElementSource/*TODO: Propagate operation info*/);
        (rslst_1,b2) = inlineReinitStmts(rslst,fns,{},false);
        true = b1 or b2;
      then
        (BackendDAE.WHEN_CLAUSE(e_1,rslst_1,io),true);
    else (inWhenClause,false);
  end matchcontinue;
end inlineWhenClause;

protected function inlineReinitStmts
"inlines function calls in reinit statements"
  input list<BackendDAE.WhenOperator> inStmts;
  input Inline.Functiontuple fns;
  input list<BackendDAE.WhenOperator> iAcc;
  input Boolean iInlined;
  output list<BackendDAE.WhenOperator> outStmts;
  output Boolean oInlined;
algorithm
  (outStmts,oInlined) := match (inStmts,fns,iAcc,iInlined)
    local
      BackendDAE.WhenOperator re;
      list<BackendDAE.WhenOperator> rest,stmts;
      Boolean b;

    case ({},_,_,_) then (listReverse(iAcc),iInlined);
    case (re::rest,_,_,_)
      equation
        (re,b) = inlineReinitStmt(re,fns);
        (stmts,b) = inlineReinitStmts(rest,fns,re::iAcc,b or iInlined);
      then
        (stmts,b);
  end match;
end inlineReinitStmts;

protected function inlineReinitStmt
"inlines function calls in a reinit statement"
  input BackendDAE.WhenOperator inReinitStatement;
  input Inline.Functiontuple inElementList;
  output BackendDAE.WhenOperator outReinitStatement;
  output Boolean inlined;
algorithm
  (outReinitStatement,inlined) := matchcontinue(inReinitStatement,inElementList)
    local
      Inline.Functiontuple fns;
      DAE.ComponentRef cref;
      DAE.Exp e,e_1;
      BackendDAE.WhenOperator rs;
      DAE.ElementSource source;
      list<DAE.Statement> assrtLst;
    case (BackendDAE.REINIT(cref,e,source),fns)
      equation
        (e_1,source,true,_) = Inline.inlineExp(e,fns,source);
      then
        (BackendDAE.REINIT(cref,e_1,source),true);
    case (rs,_) then (rs,false);
  end matchcontinue;
end inlineReinitStmt;

annotation(__OpenModelica_Interface="backend");
end BackendInline;
