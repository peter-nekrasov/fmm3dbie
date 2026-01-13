function val = v2v_adap(xt,yt)

val = integral2(@(r,th) r.*v2v_integral(r,th,xt,yt), 0, 1, 0, 2*pi);


end


function val = v2v_integral(r,th,xt,yt)
x = r.*cos(th);
y = r.*sin(th);

rs = sqrt((xt-x).^2+(yt-y).^2);
val = rs.^2.*log(rs).*test_fun(x,y)/8/pi;

end


function f = test_fun(x,y)

f = exp( - 10*x.^2 - 10*y.^2 );


end