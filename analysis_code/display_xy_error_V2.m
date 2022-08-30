function h=display_xy_error_V2(x,y, xe, ye, bar_col)


for i=1:length(x)
   
    h=plot( [x(i),x(i)], [y(i)-ye(i), y(i)+ye(i)], 'color', bar_col);
    if ~isempty(xe)
    h=plot( [x(i)-xe(i),x(i)+xe(i)], [y(i), y(i)], 'color', bar_col);
    end
     
end


end