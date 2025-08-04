function K23 = pairing_function(qt_tmp)
K1 = 0.5*(qt_tmp(1)+qt_tmp(2))*(qt_tmp(1)+qt_tmp(2)+1) + qt_tmp(2);
K12 = 0.5*(K1+qt_tmp(3))*(K1+qt_tmp(3)+1) + qt_tmp(3);
K23 = 0.5*(K12+qt_tmp(3))*(K12+qt_tmp(4)+1) + qt_tmp(4);
end