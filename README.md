# kmat

Kmat is a Ruby gem for matrix operations. Kmat uses BLAS/LAPACK as back-end.

## Requirements

BLAS/LAPACK libraries are needed to build kmat. MKL can also be used.
You need to modify extconf.rb to use the other BLAS/LAPACK compatible libraries (pull requests are welcome).

## Installation

Kmat is available on RubyGems.org:

    $ gem install kmat

or in your Gemfile:

```ruby
gem 'kmat'
```

## Usage

While the gem name is kmat, the top module of this gem is the class `Mat`. Most of the functions of kmat are defined under `Mat` (some monkey patches (e.g. `Random#randn`) are applied to built-in modules/classes).

`Mat` is a class of matricies. It has matrix operations as methods. Unlike numpy in Python, `Mat` cannot be a tensor with order other than 2. `Mat` with shape `(n, 1)` or `(1, n)` behaves as a vector (not be distinguished from each other) in some methods (e.g. `Mat#iprod`).

Matlab style expected return value comments appear in the following sample codes, but they are only in documents (kmat has no way to generate these strings from Mat instances nor way to generate Mat instances from these strings).

### Construction of matricies

```ruby
require 'kmat'

Mat[[1, 2], [3, 4]]       #=> [1, 2; 3, 4]
Mat.ones(2, 2)            #=> [1, 1; 1, 1]
Mat.eye(2, 2)             #=> [1, 0; 0, 1]
Mat.randn(2, 2)           #=> [-0.11, 0.25; 0.83, -0.03] (vary every call)
Mat.range(3)              #=> [0; 1; 2]
Mat.new(2, 3){|i, j| i+j} #=> [0, 1, 2; 1, 2, 3]
```

### Arithmetic operations
```ruby
require 'kmat'

a = Mat.ones; b = Mat[[2, 0], [0, 1]]
a+b        #=> [3, 1; 1, 2]
a-b        #=> [-1, 1; 1, 0]
a.mprod(b) #=> [2, 1; 2, 1] (matrix product)
a.e_mul(b) #=> [2, 0; 0, 1] (element-wise product)
b.under(a) #=> [0.5, 0.5; 1, 1] (like b\a in Matlab)
a.over(b)  #=> [0.5, 1; 0.5, 1] (like a/b in Matlab)
a*b        #=> ArgumentError (to avoid confusion of Mat#mprod vs Mat#e_mul)
using Mat::MatrixProductOperator
[a*b, a/b] #=> [ [2, 1; 2, 1], [0.5, 1; 0.5, 1] ] (refinements are available)
a.sin      # Most of mathematical functions defined in math.h are available as element-wise operations
```

### Destructive operation
```ruby
require 'kmat'

a = Mat.ones; b = Mat[[2, 0], [0, 1]]
a+b; a     #=> [1, 1; 1, 1]
a.add!(b)
a          #=> [3, 1; 1, 2]
c = Mat.new(2, 2)
c.mprod!(a, b)
c          #=> [6, 1; 2, 2]
a.sub!(b); a.e_mul!(b); b.e_div!(a); c.under!(a, b)
```

### Numpy-like broadcasting
```ruby
require 'kmat'

a = Mat.eye(2, 2); b = Mat.range(2)
a+b                    #=> [1, 1; 2, 2]
a.add!(a.broadcast(b)) # For destructive operation, you need to fit shape using Mat#broadcast
a.s_add!(1.5)          # Or use Mat#s_xxx! for element-wise opertion with broadcasted scalar value
0.5 * a                # Numeric#* can be used as scalar multiplication (others like Numeric+Mat are not available)
```

### Element or submatrix accessing
```ruby
rquire 'kmat'

a = Mat[[1, 2], [3, 4]]
a[1, 0]                  #=> 3.0
a[0, -1]                 #=> 2.0 (negative indexing is available)
a[nil, 1]                #=> [2; 4] (nil works like `:' in numpy)
a[[1, 0], 0..1]          #=> [3, 4; 1, 2]
a[a.gt(2)]               #=> [3; 4] (boolean indexing returns a vector)
vidx = Mat.new(1, 2, :object){ |i, j| [j, i] }
a[vidx]                  #=> [1, 3] (=[a[0, 0], a[1, 0]])
b = a[0, nil]            # a submatrix is like a `view' of numby
a[0, 0] = 5
b                        #=> [5, 2]
b[0, 1] = -1             # but it's writable and may change the `supermatrix'
a                        #=> [5, -1; 3, 4]
b.freeze                 # use Object#freeze to avoid it
a[0, 0] = -3             # but, freezing submatrix do not freeze the `supermatrix'
b                        #=> [-3, -1]
a.deep_freeze            # use `supermatrix'.deep_freeze to freeze all the related matricies
b[0, 0] = 7              #=> FrozenError
c = a[0, nil]            # submatrix made from frozen matrix is pre-frozen
c[0, 0] = 9              #=> FrozenError
a = Mat[[1, 2], [3, 4]]
a[nil, 0] = Mat.range(2) # multi-entry substitution is also available
a                        #=> [0, 2; 1, 4]
a.map(&:-@)              #=> [0, -2; -1, -4]
a.map_with_index!{|e, i, j| e+j}
a.diag                   #=> [0; 4] (returns diagonal elements as a vector)
a[]                      #=> [0, 2; 1, 4] (with no argument is equivalent to a[nil, nil])
Mat.range(3)[2]          #=> 2.0 (with single integer or single integer array is available for vectors)
```

### Value types
Mat in kmat has 5 value types.
#### Float
The default value type is `:float`. Values are `double` in C language and it is compatible with `Float` in Ruby. Most of linear argebraic operations are available only on float matricies.

#### Complex
`:complex` is available to deal with complex floats. Some operations are defined but the number of available operations is limitted.

#### Int
`:int` can be used as row or column index array. Values are `int` in C language, so it is not useful for matrix operations with large integer elements. We recomend to use `:object` with `Integer` for such usage.

#### Bool
`:bool` can be used as boolean indexing. Logical operations are available. In some operations, elements in boolean matricies behaves as a finite filed with order 2 (`+` is exclusive disjunction and `*` is logical conjunction).

#### Object
`:object` matricies can contain arbitrary Ruby objects. Operation behaviors are dependent on methods defined for the objects. For example, `Mat#solve` works if `K#==`, `K#quo`, `K#*` and `K#-` are defined appropriately, where `K` is a class of the elements.
`:object` matricies with 2-length `Array`s are used as indecies for other `Mat`s.

```ruby
rquire 'kmat'

Mat.new(2, 2, :float){|i, j| i+j}
Mat.new(2, 2, :complex){|i, j| Complex.rect(i, j)}
Mat.new(2, 2, :int){|i, j| i-j}
Mat.new(2, 2, :bool){|i, j| i==j}
Mat.new(2, 2, :object){|i, j| Rational(i, j) }
Mat.new(1, 1).vtype   #=> :float (Mat#vtype returns its value type as a Symbol above)
```

### Sort, stacking and logical operations
```ruby
require 'kmat'

a = Mat[[3, 2, 1], [5, -3, 7]]
a.sort(:row)  #=> [1, 2, 3; -3, 5, 7]
a.rsort(:col) #=> [5, 2, 7; 3, -3, 1]
a.flip(:both) #=> [7, -3, 5; 1, 2, 3]
a.t           #=> [3, 5; 2, -3; 1, 7] (transpose)
b = Mat.range(2)
Mat.blocks([[a, b], [b, a]]) #=> [3, 2, 1, 0; 5, -3, 7, 1; 0, 3, 2, 1; 1, 5, -3, 7]
Mat.vstack(a, a); Mat.hstack(a, b)
a.gt(b)       #=> [true, true, true; true, false, true] (numpy-like broadcasting)
a.eq(b); a.ne(b); a.ge(b); a.lt(b); a.le(b)
a.max         #=> 7
a.maximum(b)  #=> [3, 2, 1; 5, 1, 7]
Mat.maximum(a, b.repmat(1, 3), Mat.randn(2, 3))
              #=> Mat.maximum can recieve arbitrary number of arguments but do not broadcast automatically
```

### Linear algebraic operations
Most of them are useing BLAS/LAPACK as a back-end.
```ruby
require 'kmat'

a, b = Mat.randn(3, 3), Mat.rand(3, 1)
a.ls(b)          # ls means least_squares, the return x mimimize a.mprod(x)-b
a.evd            # eigendecomposition
a.svd            # singular value decomposition
a.lup            # LUP decomposition
a.det            # determinant
a.qr             # QR decomposition
Mat.rand_orth(3) # Random orthogonal matrix
a.tr             # trace
```

### Copying elements
```ruby
require 'kmat'

a, b = Mat.randn(2, 2), Mat.randn(2, 2)
a.dup           # dup/clone are like Array#dup/clone
a.copy_from(b)  # it is equivalent to a[] = b, value type is dependent on the reciever
a.replace(b[])  # Mat#replace copies submatrix relation and value type (Mat#copy_from do not)
a.fill(1.5)     #=> [1.5, 1.5; 1.5, 1.5] (destructive)
Marshal.load(Marshal.dump(a))
                # Marshal.dump/load are available
```

### Random
```ruby
require 'kmat'

Mat.randn(2, 2, random: Random.new) # In kmat, methods using random value, the generator can be specified by keyword `random'
$MatRandom = Random.new             # Or substitute into $MatRandom (default value of $MatRandom is Random::DEFAULT)
Random.new.randn                    # Random#randn returns a random value following N(0, 1)
randn()                             # Kernel#randn is equivalent to Random::DEFAULT.randn
```


## Contributing

Bug reports and pull requests are welcome on GitHub at https://github.com/cycloawaodorin/kmat.

## License
Kmat is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version3 of the License, or (at your option) any later version.

Kmat is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with the gem; see the file LICENSE.md. If not, see the website of [GPL](https://www.gnu.org/licenses/).

## `/ext/lapack_headers/*.h`
The header files under `/ext/lapack_headers/` directory are modified or copied version of an Azalea Clive's product or Intel Corp's products. The original versions are distributed at the following.

`blas.h` is distributed at [BLASの簡単な使い方](http://azalea.s35.xrea.com/blas/sample1.html).

`lapacke.h`, `lapacke_config.h`, `lapacke_mangling.h` and `lapacke_utils.h` are distributed at the website of [LAPACKE](http://www.netlib.org/lapack/#_standard_c_language_apis_for_lapack) under modified BSD license.
