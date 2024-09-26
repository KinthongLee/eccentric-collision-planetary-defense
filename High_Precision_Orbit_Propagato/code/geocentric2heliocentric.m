function [r_heliocentric, r_earth] = geocentric2heliocentric(Eph_eci,Mjd_UTC)
global eopdata
    r_heliocentric = zeros(size(Eph_eci,1),3);
    r_earth = zeros(size(Eph_eci,1),3);
    for i = 1 : size(Eph_eci,1)
        t = Eph_eci(i,1);
        MJD_UTC = Mjd_UTC + t/86400;
        [x_pole,y_pole,UT1_UTC,LOD,dpsi,deps,dx_pole,dy_pole,TAI_UTC] = IERS(eopdata,MJD_UTC,'l');
        [UT1_TAI,UTC_GPS,UT1_GPS,TT_UTC,GPS_UTC] = timediff(UT1_UTC,TAI_UTC);
        MJD_UT1 = MJD_UTC + UT1_UTC/86400;
        MJD_TT  = MJD_UTC + TT_UTC/86400;
        MJD_TDB = Mjday_TDB(MJD_TT);
        [~,~,~,~,~,~,~, ...
         ~,~,~,r_Sun,~] = JPL_Eph_DE440(MJD_TDB);



        r_geocentric = Eph_eci(i,2:4);
        for j = 1 : 3
            r_heliocentric(i,j) = r_geocentric(j) - r_Sun(j);
            r_earth(i,j) = -r_Sun(j);
        end
        
    end
    
end