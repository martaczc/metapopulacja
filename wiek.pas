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
 deltaas=0.01; //?
 deltaar=-0.01; //?
 czas=200;

Type population=array[1..czas+1] of real;

Var pop0, pop1: population;
    t, i, j, k, koh: longint;
    as, ar, bs, br, ps,pr,N, age:real;
    avgAge,sdAge,avgN,sdN: array[0..4] of real;
    winter: boolean;
    INFO: text;

begin
assign(INFO,'wiek.txt');
rewrite(INFO);
append(INFO);
writeln(INFO,'as ','ar ','bs ','br ','t ','avgN[t] ','sdN[t] ','avgAge[t] ','sdAge[t] ');
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
  writeln(i,' ',j ,' pr(Ne)=',1/(1+exp(-(ar*N+br))):7:3,' ps(Ne)=',1/(1+exp(-(as*N+bs))):7:3);
  for t:=0 to 4 do
   begin
   avgN[t]:=0;
   sdN[t]:=0;
   avgAge[t]:=0;
   sdAge[t]:=0;
   end;
  for t:=1 to czas do
   begin
   if ((t mod 5 = 4) or (t mod 5 = 0)) then winter:=true else winter:=false;
   age:=0;
   pr:=1/(1+exp(-(ar*N+br)));
   ps:=1/(1+exp(-(as*N+bs)));
   if (not winter) then pop1[1]:=pr*L*N/2 else pop1[1]:=0;
   N:=pop1[1];
   for k:=1 to koh do
    begin
    pop1[k+1]:=pop0[k]*ps;
    N:=N+pop1[k+1];
    age:=age+k*pop1[k+1];
    end;
   koh:=koh+1;
   if N<>0 then age:=age/N else age:=0;
   if t>100 then 
    begin
    avgN[t mod 5]:=avgN[t mod 5]+N;
    sdN[t mod 5]:=sdN[t mod 5]+N*N;
    avgAge[t mod 5]:=avgAge[t mod 5]+age;
    sdAge[t mod 5]:=sdAge[t mod 5]+age*age;
    end;
   for k:=1 to koh do pop0[k]:=pop1[k];
   //pop0:=pop1;
   gotoxy(1,whereY);
   write(i,' ',j,' ',t)
   end;
  for t:=0 to 4 do
   begin
   avgN[t]:=avgN[t]/(czas-100);
   sdN[t]:=sqrt(sdN[t]/(czas-100) - avgN[t]*avgN[t]);
   avgAge[t]:=avgAge[t]/(czas-100);
   sdAge[t]:=sqrt(sdAge[t]/(czas-100) - avgAge[t]*avgage[t]);
   append(INFO);
   writeln(INFO,as:7:2,' ',ar:7:2,' ',bs:7:5,' ',br:7:5,' ',t,' ',avgN[t]:7:5,' ',sdN[t]:7:5,' ',avgAge[t]:7:5,' ',sdAge[t]:7:5);
   close(INFO);
   end;
  end;
 end;
end.
