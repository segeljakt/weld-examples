# Weld IR Examples

Below are algorithms implemented in [Weld](https://github.com/weld-project/weld).

## Matrix Multiplication

```
let n = 2L;
let p = 2L;
let m = 3L;
# nxm
let A = [1,2,3,
         4,5,6];
# mxp
let B = [1,2,
         4,5,
         7,8];
# nxp
let C = # AxB
  for(rangeiter(0L, n, 1L), appender[i32], |C,x,i|
    for(rangeiter(0L, p, 1L), C, |C,x,j|
      let sum =
        for(rangeiter(0L, m, 1L), merger[i32,+], |sum,x,k|
          let Aik = lookup(A, i*n+k);
          let Bkj = lookup(B, k*m+j);
          merge(sum, Aik * Bkj)
        );
      merge(C, result(sum))
    )
  );
result(C)
```

## PageRank 

```
let src = [0L,0L,1L,2L];
let dst = [1L,2L,2L,0L];

let out_edges = zip(src, dst);
let in_edges  = zip(dst, src);

# out_edges   in_edges
# [{0,1},     [{1,0},
#  {0,2},      {2,0},
#  {1,2},      {2,1},
#  {2,0}]      {0,2}]

let in_nbrs = groupmerger[i64,i64];
let in_nbrs = for(in_edges, in_nbrs, |b,i,x| merge(b, x));
let in_nbrs = result(in_nbrs);

# ({0,[2]},
#  {1,[0]},
#  {2,[0,1]})

let fan_outs = dictmerger[i64,i64,+];
let fan_outs = for(out_edges, fan_outs, |b,i,x| merge(b, {x.$0, 1L}));
let fan_outs = result(fan_outs);

# ({0,2},
#  {1,1},
#  {2,1})

let n = merger[i64,max];
let n = for(src, n, |b,i,x| merge(b, x));
let n = for(dst, n, |b,i,x| merge(b, x));
let n = result(n)+1L;

# n = 3

let initial_rank = 1.0/f64(n);
let initial_ranks = appender[f64];
let initial_ranks = for(rangeiter(0L, n, 1L), initial_ranks, |b,i,x| merge(b, initial_rank));
let initial_ranks = result(initial_ranks);

# [0.33,
#  0.33,
#  0.33]

let start = {initial_ranks, 0};

let teleport       = 0.1;
let tolerance      = 0.0001;
let max_iterations = 20;

iterate(start, |iterator|
  let old_ranks = iterator.$0;
  let iteration = iterator.$1;

  # PageRank
  let ranks = result(for(rangeiter(0L, n, 1L), appender[f64], |ranks,i,node|
    let sum = result(for(lookup(in_nbrs, node), merger[f64,+], |sum,j,in_nbr|
      let fan_out = lookup(fan_outs, in_nbr);
      let old_rank = lookup(old_ranks, in_nbr);
      merge(sum, old_rank / f64(fan_out))
    ));
    merge(ranks, teleport/f64(n) + (1.0 - teleport) * sum)
  ));

  # Condition
  let max_delta = result(for(zip(ranks,old_ranks), merger[f64,max], |b,i,x|
    let delta = x.$0-x.$1;
    let abs = select(delta < 0.0, -1.0*delta, delta);
    merge(b, abs)
  ));
  let cond1 = select(iteration < max_iterations, true, false);
  let cond2 = select(max_delta < tolerance, false, true);
  {{ranks, iteration+1}, cond1 && cond2}
)

# PageRank:
#
# [0.388,
#  0.215,
#  0.397]
```
