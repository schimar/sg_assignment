{
ref = length($4);
alt = length($5); 
if (ref == 1 && alt == 1) snp ++ ; 
else if (ref > 1 && alt == 1) del ++ ;
else ins ++
}
END{
print "SNP: ", snp, "\nDEL: ", del, "\nINS: ", ins;
}
