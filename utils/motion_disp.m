function motion_disp(errs, LK, id)
if id >0
    errs = errs([1:id-1,id+1:end]);
    LK  = LK([1:id-1,id+1:end]);
end

err2obj = mean(errs(LK==2));
mederr2obj = median(errs(LK==2));

err3obj = mean(errs(LK==3));
mederr3obj = median(errs(LK==3));

errobj =   mean(errs);
mederrobj = median(errs);

disp([err2obj, mederr2obj;
    err3obj, mederr3obj;
    errobj, mederrobj ] );

end