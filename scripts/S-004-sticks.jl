gr()

sticks = [0 1 0 1 0 1 0
          1 0 1 0 1 0 1
          0 1 0 1 0 1 0
          1 0 1 0 1 0 1
          0 1 0 1 0 1 0
          1 0 1 0 1 0 1
          0 1 0 1 0 1 0];
#
test_sticks = [0 1 0 1 0 1 0
               1 0 1 0 1 0 1
               0 1 0 0 0 1 0
               1 0 1 0 1 0 1
               0 1 0 1 0 1 0
               1 0 1 0 1 0 1
               0 1 0 0 0 1 0]

sum(sticks)
sum(test_sticks)

using Plots

#            _|              _|
#  _|_|_|    _|    _|_|    _|_|_|_|
#  _|    _|  _|  _|    _|    _|
#  _|    _|  _|  _|    _|    _|
#  _|_|_|    _|    _|_|        _|_|
#  _|
#  _|


function plot_sticks(sticks)
  plot(size=(500,450), xlim=(0.5,4.5), ylim=(0.5,4.5))
  Δ = 0.1
  for i in 1:7, j=1:7
    if sticks[8-i,j]==1
      ii = div(i+1,2)
      jj = div(j+1,2)
      if isodd(i)
        plot!([jj+Δ,jj+1-Δ], [ii,ii], lab="", lw=6, c=:black)
        scatter!([jj+1-Δ], [ii], m=(7,:o,:red), lab="")
      end
      if iseven(i)
        plot!([jj,jj], [ii+Δ,ii+1-Δ],lab="", lw=6, c=:black)
        scatter!([jj], [ii+Δ], m=(7,:o,:red), lab="")
      end
    end
  end
  plot!(xaxis=false, yaxis=false)
end

plot_sticks(sticks)
plot_sticks(test_sticks)


#      _|_|                                  _|      _|
#    _|      _|    _|  _|_|_|      _|_|_|  _|_|_|_|        _|_|    _|_|_|      _|_|_|
#  _|_|_|_|  _|    _|  _|    _|  _|          _|      _|  _|    _|  _|    _|  _|_|
#    _|      _|    _|  _|    _|  _|          _|      _|  _|    _|  _|    _|      _|_|
#    _|        _|_|_|  _|    _|    _|_|_|      _|_|  _|    _|_|    _|    _|  _|_|_|


function count_small_squares(sticks)
  counter = 0
  for i in [2,4,6], j in [2,4,6]
    if sticks[i,j+1]+sticks[i+1,j]+sticks[i,j-1]+sticks[i-1,j] == 4
      counter += 1
    end
  end
  return counter
end

function count_medium_squares(sticks)
  counter = 0
  for i in [3,5], j in [3,5]
    if sticks[i-2,j+1]+sticks[i-2,j-1]+
        sticks[i-1,j-2]+sticks[i+1,j-2]+
            sticks[i+2,j-1]+sticks[i+2,j+1]+
                sticks[i+1,j+2]+sticks[i-1,j+2] == 8
      counter += 1
    end
  end
  return counter
end

function count_big_squares(sticks)
  counter = 0
  if sum(sticks[1,:])+ sum(sticks[end,:])+ sum(sticks[:,1])+ sum(sticks[:,end]) == 12
    counter += 1
  end
  return counter
end

#
#    _|                            _|
#  _|_|_|_|    _|_|      _|_|_|  _|_|_|_|    _|_|_|
#    _|      _|_|_|_|  _|_|        _|      _|_|
#    _|      _|            _|_|    _|          _|_|
#      _|_|    _|_|_|  _|_|_|        _|_|  _|_|_|


count_all_squares(sticks) =
  count_small_squares(sticks) +
  count_medium_squares(sticks) +
  count_big_squares(sticks)

sticks
plot_sticks(sticks)
count_small_squares(sticks)
count_medium_squares(sticks)
count_big_squares(sticks)
count_all_squares(sticks)
#
test_sticks
plot_sticks(test_sticks)
count_small_squares(test_sticks)
count_medium_squares(test_sticks)
count_big_squares(test_sticks)
count_all_squares(test_sticks)

test3_sticks = [0 1 0 1 0 0 0
                1 0 1 0 1 0 0
                0 1 0 0 0 0 0
                1 0 0 0 1 0 0
                0 1 0 1 0 0 0
                1 0 1 0 0 0 0
                0 1 0 0 0 0 0]

#
plot_sticks(test3_sticks)
count_small_squares(test3_sticks)
count_medium_squares(test3_sticks)
count_big_squares(test3_sticks)
count_all_squares(test3_sticks)


#                                          _|    _|
#  _|  _|_|    _|_|      _|_|_|  _|    _|  _|  _|_|_|_|    _|_|_|
#  _|_|      _|_|_|_|  _|_|      _|    _|  _|    _|      _|_|
#  _|        _|            _|_|  _|    _|  _|    _|          _|_|
#  _|          _|_|_|  _|_|_|      _|_|_|  _|      _|_|  _|_|_|


function remove_one_stick(sticks)
  new_sticks = copy(sticks)
  #
  while true
    i,j = rand(1:7), rand(1:7)
    if sticks[i,j]==1
      new_sticks[i,j] = 0
      return new_sticks
    end
  end
end

function remove_four(sticks)
  v = remove_one_stick(sticks)
  v = remove_one_stick(v)
  v = remove_one_stick(v)
  v = remove_one_stick(v)
  return v
end
#

plot_sticks(sticks)

samples = [remove_four(sticks) for i = 1:1_000_000]

Nsquares = count_all_squares.(samples)

histogram(Nsquares, bins=0.5:14.5)

fsamples = samples[Nsquares .== 2]

bb = 1
begin
  # fsamples[bb]
  plot_sticks(fsamples[bb])
  bb = bb+1
  plot!()
end
