% use module(library(lists))

nat(0).
nat(s(X)) :- nat(X).

plus(0, X, X).
plus(s(X), Y, s(Z)) :- plus(X, Y, Z).

minus(X, Y, Z) :- plus(Y, Z, X).

fib(0, s(0)).
fib(s(0), s(0)).
fib(s(s(K)), N) :- fib(s(K), A), fib(K, B), plus(A, B, N).

lista([]).
lista([_ | L]) :- lista(L).

pierwszy(E, [E | _]).

ostatni(E, [E]).
ostatni(E, [_ | L]) :- ostatni(E, L).

element(E, [E | _]).
element(E, [_ | L]) :- element(E, L).


scal([], L, L).
scal([A | L1], L2, [A | L3]) :- scal(L1, L2, L3).

podziel([], [], []).
podziel([A], [A], []).
podziel([A, B | L], [A | N], [B | P]) :- podziel(L, N, P).

rowna([], []).
rowna([A | L1], [A | L2]) :- rowna(L1, L2).

prefiks([], _).
prefiks([A | L1], [A | L2]) :- prefiks(L1, L2).

podlista(L1, L2) :- prefiks(L1, L2).
podlista(L1, [_ | L2]) :- prefiks(L1, L2).

podciag([], _).
podciag([A | L1], [A | L2]) :- podciag(L1, L2).
podciag(L1, [_ | L2]) :- podciag(L1, L2).

suma([], 0).
suma([A | L], A + N) :- suma(L, N).

% dlugosc0(+L, ?N)
dlugosc0([], 0).
dlugosc0([_ | L], N + 1) :- dlugosc0(L, N).

% dlugosc1(+L, +N)
dlugosc1([], 0).
dlugosc1([_ | L], N) :- N1 is N - 1, dlugosc1(L, N1).


niewieksza([], _).
niewieksza([A | L], M) :- niewieksza(L, M), M =< A.

min([A], A).
min(L, M) :- member(M, L), niewieksza(L, M).

odwroc(L,R) :- odwroc(L, [], R).
odwroc([], M, M).
odwroc([L1|LS], M, R) :- odwroc(LS, [L1|M], R).

palindrom(L) :- palindrom(L, [], []).
palindrom([], A, A).
palindrom([E | L], A, B) :- append(B, [E], B1), palindrom(L, [E | A], B1).

maa([]).
maa([a | L]) :- maa(L).

mab([]).
mab([b | L]) :- mab(L).

slowo(L) :- prefiks(P, L), append(P, R, L), maa(P), mab(R).

flagaPolska(L, F) :- flagaPolska(L, [], F).
flagaPolska([], [], []).
flagaPolska([b | L], C, [b | F]) :- flagaPolska(L, C, F).
flagaPolska([c | L], C, F) :- flagaPolska(L, [c | C], F).
flagaPolska([], [c | C], [c | F]) :- flagaPolska([], C, F).

mniejsze([], _, []).
mniejsze([A | L], E, [A | M]) :- A < E, mniejsze(L, E, M).
mniejsze([A | L], E, M) :- A >= E, mniejsze(L, E, M).

wiekszeR([], _, []).
wiekszeR([A | L], E, [A | M]) :- A >= E, wiekszeR(L, E, M).
wiekszeR([A | L], E, M) :- A < E, wiekszeR(L, E, M).

quickSort([], []).
quickSort([E | L], S) :- mniejsze(L, E, M), wiekszeR(L, E, W), quickSort(M, MS), quickSort(W, WS), append(MS, [E | WS], S).

% drzewo([]).
% drzewo((_, L, R)) :- drzewo(L), drzewo(R).
% drzewo((_, [], [])).

% insertBST([], E, [E, [], []]).
% insertBST([E, L, R], N, [E, L1, R]) :- N =< E, insertBST(L, N, L1).
% insertBST([E, L, R], N, [E, L, R1]) :- N > E, insertBST(R, N, R1).

% wypiszBST([], []).
% wypiszBST([A | L], [E, LT, RT]) :- A = E, append(L1, L2, L), wypiszBST(L1, LT), wypiszBST(L2, RT).


% stworzBST([], []).
% stworzBST([A | L], D) :- stworzBST(L, D1), insertBST(D1, A, D).

liscie([], []).
liscie([E, [], []], [E]).
liscie([_, L, R], X) :- liscie(L, LL), liscie(R, RL), append(LL, RL, X).



edge(a, b).
edge(b, c).
edge(c, d).

connect(A, B) :- edge(A, B).
connect(A, C) :- edge(A, B), connect(B, C).

path(A, A, [A]).
path(A, B, [A, P | L]) :- edge(A, P), path(P, B, [P |L]).

edgeC(A, B) :- edge(A, B).
edgeC(A, B) :- edge(B, A).

pathC(A, B, P) :- pathC(A, B, [], P).
pathC(A, B, V, [A, B]) :- edge(A, B), \+member(A, V).
pathC(A, B, V, [A, X | Z]) :-
  edge(A, X),
  \+member(A, V),
  pathC(X, B, [A | V], [X | Z]).

zawieraEl([], _).
zawieraEl([A | L], R) :- member(A, R), zawieraEl(L, R).

euler([]).
euler([_]).
euler(P) :- pathC(_, _, V, _), zawieraEl(P, V), zawieraEl(V, P).

dodaj(Elem, Lista) :- var(Lista), !, Lista = [Elem | _ ].
dodaj(Elem, [ _ | Lista]) :- dodaj(Elem, Lista).


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% tree(empty).
% tree(node(_, L, P)) :- tree(L), tree(P).

% wypiszBST([], []).
% wypiszBST([A | L], [E, LT, RT]) :- A = E, append(L1, L2, L), wypiszBST(L1, LT), wypiszBST(L2, RT).


insertBST(empty, E, node(E, empty, empty)).
insertBST(node(E, L, R), N, node(E, L1, R)) :- N @< E, insertBST(L, N, L1).
insertBST(node(E, L, R), N, node(E, L, R1)) :- N @> E, insertBST(R, N, R1).
insertBST(node(E, L, R), E, node(E, L, R)).

makeLetterBST([], empty).
makeLetterBST([fp(_, C, _) | L], T) :- makeLetterBST(L, T1), insertBST(T1, C, T).

makeStateBST([], empty).
makeStateBST([fp(S1, _, S2) | L], T) :- makeStateBST(L, T1), insertBST(T1, S1, T2), insertBST(T2, S2, T).

% TODO styl, długości linii, może bez mieszania polski/eng

treeSize(empty, 0).
treeSize(node(_, L, R), 1 + LSize + RSize) :- treeSize(L, LSize), treeSize(R, RSize).

insertStateRelationsTree(empty, fp(S1, C, S2), node((S1, T), empty, empty)) :- insertRelationTree(empty, p(C, S2), T).
insertStateRelationsTree(node((S, T), L, R), fp(S1, C, S2), node((S, T), L1, R)) :- S1 @< S, insertStateRelationsTree(L, fp(S1, C, S2), L1).
insertStateRelationsTree(node((S, T), L, R), fp(S1, C, S2), node((S, T), L, R1)) :- S1 @> S, insertStateRelationsTree(R, fp(S1, C, S2), R1).
insertStateRelationsTree(node((S, T), L, R), fp(S, C, S2), node((S, T1), L, R)) :- insertRelationTree(T, p(C, S2), T1).

insertRelationTree(empty, p(C, S), node(p(C, S), empty, empty)).
insertRelationTree(node(p(C1, S1), L, R), p(C2, S2), node(p(C1, S1), L1, R)) :- C2 @< C1, insertRelationTree(L, p(C2, S2), L1).
insertRelationTree(node(p(C1, S1), L, R), p(C2, S2), node(p(C1, S1), L, R1)) :- C2 @> C1, insertRelationTree(R, p(C2, S2), R1).


makeStateRelationsTree([], empty).
makeStateRelationsTree([E | L], T) :- makeStateRelationsTree(L, T1), insertStateRelationsTree(T1, E, T).


% correct(dfa(FP, Start, [Final | FinalList]), _) :- 
%     makeStateRelationsTree(FP, StateTree), 
%     stateInStateTree(Start, StateTree), 
%     allStatesInStateTree([Final | FinalList], StateTree),


% start in states
% all finals in states
% each state is connected to all other states with all letters ---- find state in statetree, check size = lettercount
















