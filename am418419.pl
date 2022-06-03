%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adrian Matwiejuk
% TODO zwięzły opis przyjętej wewnętrznej reprezentacji automatu.


% TODO styl, długości linii, może bez mieszania polski/eng

insertBST(empty, E, node(E, empty, empty)).
insertBST(node(E, L, R), N, node(E, L1, R)) :- N @< E, insertBST(L, N, L1).
insertBST(node(E, L, R), N, node(E, L, R1)) :- N @> E, insertBST(R, N, R1).
insertBST(node(E, L, R), E, node(E, L, R)).

makeLetterBST([], empty).
makeLetterBST([fp(_, C, _) | L], T) :- 
    makeLetterBST(L, T1), 
    insertBST(T1, C, T).

makeStateBSTFromFP([], empty).
makeStateBSTFromFP([fp(S1, _, S2) | L], T) :- 
    makeStateBSTFromFP(L, T1), 
    insertBST(T1, S1, T2), 
    insertBST(T2, S2, T).

makeStateBSTFromList([], empty).
makeStateBSTFromList([S | L], T) :-
    makeStateBSTFromList(L, T1),
    insertBST(T1, S, T).

treeSize(empty, 0).
treeSize(node(_, L, R), 1 + LSize + RSize) :- 
    treeSize(L, LSize), 
    treeSize(R, RSize).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

insertStateRelationsTree(empty, fp(S1, C, S2), node(st(S1, T), empty, empty)) :- 
    insertRelationTree(empty, p(C, S2), T).
insertStateRelationsTree(node(st(S, T), L, R), fp(S1, C, S2), node(st(S, T), L1, R)) :- 
    S1 @< S, 
    insertStateRelationsTree(L, fp(S1, C, S2), L1).
insertStateRelationsTree(node(st(S, T), L, R), fp(S1, C, S2), node(st(S, T), L, R1)) :- 
    S1 @> S, 
    insertStateRelationsTree(R, fp(S1, C, S2), R1).
insertStateRelationsTree(node(st(S, T), L, R), fp(S, C, S2), node(st(S, T1), L, R)) :- 
    insertRelationTree(T, p(C, S2), T1).

insertRelationTree(empty, p(C, S), node(p(C, S), empty, empty)).
insertRelationTree(node(p(C1, S1), L, R), p(C2, S2), node(p(C1, S1), L1, R)) :- 
    C2 @< C1, 
    insertRelationTree(L, p(C2, S2), L1).
insertRelationTree(node(p(C1, S1), L, R), p(C2, S2), node(p(C1, S1), L, R1)) :- 
    C2 @> C1, 
    insertRelationTree(R, p(C2, S2), R1).

makeStateRelationsTree([], empty).
makeStateRelationsTree([E | L], T) :- 
    makeStateRelationsTree(L, T1), 
    insertStateRelationsTree(T1, E, T).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

elementInTree(E, node(E, _, _)).
elementInTree(E1, node(E2, L, _)) :-
    E1 @< E2,
    elementInTree(E1, L).
elementInTree(E1, node(E2, _, R)) :- 
    E1 @> E2, 
    elementInTree(E1, R).

allElementsInTree([], _).
allElementsInTree([E | L], T) :- 
    elementInTree(E, T),
    allElementsInTree(L, T).
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stateFullRelations(S, node(st(S, T), _, _), LetterCount1) :- treeSize(T, LetterCount2), A is LetterCount1, A is LetterCount2.
stateFullRelations(S1, node(st(S2, _), L, _), LetterCount) :-
    S1 @< S2, 
    stateFullRelations(S1, L, LetterCount).
stateFullRelations(S1, node(st(S2, _), _, R), LetterCount) :-
    S1 @> S2, 
    stateFullRelations(S1, R, LetterCount).


allStatesFullRelation(empty, _, _).
allStatesFullRelation(node(S, L, R), STR, LetterCount) :- 
    stateFullRelations(S, STR, LetterCount),
    allStatesFullRelation(L, STR, LetterCount),
    allStatesFullRelation(R, STR, LetterCount).

% correct(+Automata, -Representation)
correct(dfa(FP, Start, Finals), idfa(StateRelationsTree, Start, FinalStateSet)) :-
    makeStateBSTFromFP(FP, StateTree),
    allElementsInTree([Start | Finals], StateTree),
    makeLetterBST(FP, LetterTree),
    treeSize(LetterTree, LetterCount),
    makeStateRelationsTree(FP, StateRelationsTree),
    allStatesFullRelation(StateTree, StateRelationsTree, LetterCount),
    makeStateBSTFromList(Finals, FinalStateSet).
% start in states
% all finals in states
% each state is connected with all letters ---- find state in statetree, check size = lettercount


stateRelations(S, node(st(S, T), _, _), T).
stateRelations(S1, node(st(S2, _), L, _), T) :-
    S1 @< S2, 
    stateRelations(S1, L, T).
stateRelations(S1, node(st(S2, _), _, R), T) :-
    S1 @> S2, 
    stateRelations(S1, R, T).

nextState(C, node(p(C, S), _, _), S).
nextState(C1, node(p(C2, _), L, _), S) :-
    C1 @< C2,
    nextState(C1, L, S).
nextState(C1, node(p(C2, _), _, R), S) :-
    C1 @> C2,
    nextState(C1, R, S).

% accept(+Automata, ?Word) TODO działa dobrze tylko dla grounded
accept(A, W) :-
    correct(A, idfa(StateRelationsTree, Start, FinalStateSet)),
    acceptNext(W, Start, StateRelationsTree, FinalStateSet).

acceptNext([], S, _, FinalStateSet) :- elementInTree(S, FinalStateSet).    
acceptNext([C | W], S, STR, FinalStateSet) :-
    stateRelations(S, STR, SR),  % get the current state's relation tree
    nextState(C, SR, NS),
    acceptNext(W, NS, STR, FinalStateSet).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% empty(+Automata)
empty(A) :-
    correct(A, AR),
    emptyInternal(AR).

emptyInternal(idfa(StateRelationsTree, Start, FinalStateSet)) :-
    \+getToFinal(Start, StateRelationsTree, FinalStateSet). % if unable to get from start to final state then language is empty

getToFinal(S, STR, FSS) :- getToFinal(S, [], STR, FSS).
getToFinal(S, _, _, FinalStateSet) :- elementInTree(S, FinalStateSet).   
getToFinal(S, V, STR, FinalStateSet) :-
    stateRelations(S, STR, SR),  % get the current state's relation tree
    nextState(_, SR, NS),
    \+member(S, V),
    getToFinal(NS, [S | V], STR, FinalStateSet).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nonFinalStateList(empty, _, SL, SL).
nonFinalStateList(node(st(S, _), L, R), FinalStateSet, SL1, SL) :-
    \+elementInTree(S, FinalStateSet),
    SL2 = [S | SL1],
    nonFinalStateList(L, FinalStateSet, SL2, SL3),
    nonFinalStateList(R, FinalStateSet, SL3, SL).
nonFinalStateList(node(st(S, _), L, R), FinalStateSet, SL1, SL) :-
    elementInTree(S, FinalStateSet),
    nonFinalStateList(L, FinalStateSet, SL1, SL2),
    nonFinalStateList(R, FinalStateSet, SL2, SL).

complement(idfa(STR, Start, FinalStateSet), idfa(STR, Start, SwappedFinalStateSet)) :-
    nonFinalStateList(STR, FinalStateSet, [], SwappedFinalStateList),
    makeStateBSTFromList(SwappedFinalStateList, SwappedFinalStateSet).

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


allLettersInTree([], _).
allLettersInTree([fp(_, C, _) | L], T) :- 
    elementInTree(C, T),
    allLettersInTree(L, T).

alphabetsEqual(dfa(FP1, _, _), dfa(FP2, _, _)) :-
    makeLetterBST(FP1, LetterTree1),
    allLettersInTree(FP2, LetterTree1),
    makeLetterBST(FP2, LetterTree2),
    allLettersInTree(FP1, LetterTree2).

stateRelationsTreeToStateList(T, L) :- stateRelationsTreeToStateList(T, [], L).
stateRelationsTreeToStateList(empty, SL, SL).
stateRelationsTreeToStateList(node(st(S, _), L, R), SL, SL3) :-
    SL1 = [S | SL],
    stateRelationsTreeToStateList(L, SL1, SL2),
    stateRelationsTreeToStateList(R, SL2, SL3).

addCrossProduct(_, [], L, L).
addCrossProduct(A, [B | BL], CL, [(A, B) | CL1]) :- addCrossProduct(A, BL, CL, CL1).

makeCrossProduct(A, B, C) :- makeCrossProduct(A, B, [], C).
makeCrossProduct([], _, L, L).
makeCrossProduct([A | AL], BL, CL, CL2) :- 
    makeCrossProduct(AL, BL, CL, CL1),
    addCrossProduct(A, BL, CL1, CL2).

treeToList(T, L) :- treeToList(T, [], L).
treeToList(empty, L, L).
treeToList(node(E, LT, RT), L, L3) :-
    L1 = [E | L],
    treeToList(LT, L1, L2),
    treeToList(RT, L2, L3).

stateCrossProduct(STR1, STR2, SL) :-
    stateRelationsTreeToStateList(STR1, SL1),
    stateRelationsTreeToStateList(STR2, SL2),
    makeCrossProduct(SL1, SL2, SL).

finalStateCrossProduct(FSS1, FSS2, FSS) :-
    treeToList(FSS1, FSL1),
    treeToList(FSS2, FSL2),
    makeCrossProduct(FSL1, FSL2, FSL),
    makeStateBSTFromList(FSL, FSS).

addIntersectionRelations([], _, _, _, _, RL, RL).
addIntersectionRelations([C | Alphabet], S1, SR1, S2, SR2, RL1, RL) :-
    nextState(C, SR1, T1),
    nextState(C, SR2, T2),
    RL2 = [fp((S1, S2), C, (T1, T2)) | RL1],
    addIntersectionRelations(Alphabet, S1, SR1, S2, SR2, RL2, RL).

addAllIntersectionRelations(Alphabet, SL, STR1, STR2, FP) :- addAllIntersectionRelations(Alphabet, SL, STR1, STR2, [], FP).
addAllIntersectionRelations(_, [], _, _, RL, RL).
addAllIntersectionRelations(Alphabet, [(S1, S2) | SL], STR1, STR2, RL1, RL) :-
    stateRelations(S1, STR1, SR1),
    stateRelations(S2, STR2, SR2),
    addIntersectionRelations(Alphabet, S1, SR1, S2, SR2, RL1, RL2),
    addAllIntersectionRelations(Alphabet, SL, STR1, STR2, RL2, RL).

intersection(Alphabet, idfa(STR1, Start1, FSS1), idfa(STR2, Start2, FSS2), idfa(STR, (Start1, Start2), FSS)) :-
    finalStateCrossProduct(FSS1, FSS2, FSS),
    stateCrossProduct(STR1, STR2, SL),
    addAllIntersectionRelations(Alphabet, SL, STR1, STR2, FP),
    makeStateRelationsTree(FP, STR).

alphabet(dfa(FP, _, _), Alphabet) :-
    makeLetterBST(FP, LetterTree),
    treeToList(LetterTree, Alphabet).



% A = 0, 2, 4, 6, 8
% A'= 1, 3, 5, 7, 9
% B = 0, 4, 8
% B'= 1, 2, 3, 5, 6, 7, 9

% a3 = a5 <=>  A = B <=> A n B' = (2,6,...) = empty FALSE
% a5 = a3 <=>  B = A <=> B n A' = () = empty TRUE


% subsetEq(+Automata1, +Automata2) 
% A c B <=> (A n B') = empty
subsetEq(A, B) :-
    alphabetsEqual(A, B),
    correct(A, AR),
    correct(B, BR),
    alphabet(A, Alphabet),
    subsetInternal(Alphabet, AR, BR).

% AR c BR <=> (AR n BR') = empty
subsetInternal(Alphabet, AR, BR) :-
    complement(BR, CR),
    intersection(Alphabet, AR, CR, I),
    emptyInternal(I).

% equal(+Automata1, +Automata2)
% A = B <=> A c B and B c A
equal(A, B) :-
    alphabetsEqual(A, B),
    correct(A, AR),
    correct(B, BR),
    alphabet(A, Alphabet),
    subsetInternal(Alphabet, AR, BR),
    subsetInternal(Alphabet, BR, AR).


example(a11, dfa([fp(1,a,1),fp(1,b,2),fp(2,a,2),fp(2,b,1)], 1, [2,1])).
example(a12, dfa([fp(x,a,y),fp(x,b,x),fp(y,a,x),fp(y,b,x)], x, [x,y])).
example(a2, dfa([fp(1,a,2),fp(2,b,1),fp(1,b,3),fp(2,a,3), fp(3,b,3),fp(3,a,3)], 1, [1])).
example(a3, dfa([fp(0,a,1),fp(1,a,0)], 0, [0])).
example(a4, dfa([fp(x,a,y),fp(y,a,z),fp(z,a,x)], x, [x])).
example(a5, dfa([fp(x,a,y),fp(y,a,z),fp(z,a,zz),fp(zz,a,x)], x, [x])).
example(a6, dfa([fp(1,a,1),fp(1,b,2),fp(2,a,2),fp(2,b,1)], 1, [])).
example(a7, dfa([fp(1,a,1),fp(1,b,2),fp(2,a,2),fp(2,b,1),fp(3,b,3),fp(3,a,3)], 1, [3])).

% bad ones
example(b1, dfa([fp(1,a,1),fp(1,a,1)], 1, [])).
example(b2, dfa([fp(1,a,1),fp(1,a,2)], 1, [])).
example(b3, dfa([fp(1,a,2)], 1, [])).
example(b4, dfa([fp(1,a,1)], 2, [])).
example(b5, dfa([fp(1,a,1)], 1, [1,2])).
example(b6, dfa([], [], [])).
