{program liczacy spodziewana strukture wiekowa populacji w zaleznosci od wartosci parametrow rozrodu, przy ustalonym zageszczeniu rownowagi i prawdopodobienstwie smierci w stanie rownowagi (model deterministyczny)}

{coś jest nie tak: za mała liczebność populacji, za duże sd wieku i liczebnosci}


Uses Crt;
const
 L=6;//6 mlodych
 psNe=0.6;
 Ne=20;
 nas=100; //liczba sprawdzanych wartosci as
 nar=100; //liczba sprawdzanych wartosci ar
 as0=0;//?
 ar0=0;//?
 deltaas=0.002; //?
 deltaar=-0.002; //?
 czas=200;

Type population=array[1..(czas+1)] of real;

Var pop0, pop1: population;
    t, i, j, k, koh: longint;
    as, ar, bs, br, ps,pr,N, age:real;
    avgAge,varAge,avgN,varN: array[0..4] of real;
    winter: boolean;
    INFO: text;

begin
assign(INFO,'wiek.txt');
rewrite(INFO);
append(INFO);
writeln(INFO,'as ','ar ','bs ','br ','t ','avgN[t] ','varN[t] ','avgAge[t] ','varAge[t]');
close(INFO);
for i:=0 to nas do
 begin
 as:=as0+i*deltaas;
 bs:= ln(psNe/(1-psNe)) - as*Ne; 
 for j:=0 to nar do
  begin
  ar:=ar0+j*deltaar;
  br:= ln(2*psNe/(L-2*psNe)) - ar*Ne;
  koh:=1;
  pop0[1]:=20;
  N:=20;
  for t:=0 to 4 do
   begin
   avgN[t]:=0;
   varN[t]:=0;
   avgAge[t]:=0;
   varAge[t]:=0;
   end;
  for t:=1 to czas do
   begin
   if ((t mod 5 = 4) or (t mod 5 = 0)) then winter:=true else winter:=false;
   //winter:=false;
   age:=0;
   pr:=1/(1+exp(-(ar*N+br)));
   ps:=1/(1+exp(-(as*N+bs)));
   if (not winter) then pop1[1]:=pr*L*N/2 else pop1[1]:=0;
   for k:=1 to koh do pop1[k+1]:=pop0[k]*(1-ps);
   koh:=koh+1;
   N:=0;
   for k:=1 to koh do
    begin 
    pop0[k]:=pop1[k];
    age:=age+k*pop1[k];
    N:=N+pop1[k];
    end;
   if N<>0 then age:=age/N else age:=0;
   if t>100 then 
    begin
    avgN[t mod 5]:=avgN[t mod 5]+N;
    varN[t mod 5]:=varN[t mod 5]+N*N;
    avgAge[t mod 5]:=avgAge[t mod 5]+age;
    varAge[t mod 5]:=varAge[t mod 5]+age*age;
    end;
   gotoxy(1,whereY);
   write(i,' ',j,' ',t);
   end;
  append(INFO);
  for t:=0 to 4 do
   begin
   avgN[t]:=avgN[t]/((czas-100)/5);
   varN[t]:=varN[t]/((czas-100)/5) - avgN[t]*avgN[t];
   avgAge[t]:=avgAge[t]/((czas-100)/5);
   varAge[t]:=varAge[t]/((czas-100)/5) - avgAge[t]*avgage[t];
   writeln(INFO,as:7:5,' ',ar:7:5,' ',bs:7:5,' ',br:7:5,' ',t,' ',avgN[t]:7:5,' ',varN[t]:7:5,' ',avgAge[t]:7:5,' ',varAge[t]:7:5);
   end;
  close(INFO);
  end;
 end;
end.
