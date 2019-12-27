% for k=33 add 5mK to DOS for 8.5 channels

k=40
DDDD=DDD;
for c=k+1:110000
   DDDD(c)=DDD(c-k);
end
DDD=DDDD;
save
