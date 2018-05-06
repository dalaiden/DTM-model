function[lot_year_one_data] = corres_data_nb_years

periode_years_ans = 5;
nb_annee = 6;
lot_year_one_data = nan(periode_years_ans*nb_annee*360,2);

debut_fin_annee = nan(nb_annee,2);
debut_fin_annee(1,1) = 1;
debut_fin_annee(1,2) = nb_annee*360;

for i = 2:periode_years_ans;
    debut_fin_annee(i,1) = debut_fin_annee(i-1,2) + 1;
    debut_fin_annee(i,2) = debut_fin_annee(i-1,2) + nb_annee*360;
end 

for i = 1:length(lot_year_one_data);
    lot_year_one_data(i,1) = i;
end 

for i = 1:periode_years_ans;
    lot_year_one_data(debut_fin_annee(i,1):debut_fin_annee(i,2),2) = 1:nb_annee*360;
end