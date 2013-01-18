//Kod źródłowy programu symulującego metapopulację złożoną z subpopulacji żyjących na rozłącznych płatach powierzchni – model 19d
{dodać:
-jak zgrabnie ustawić zmienną liczbę genów?
-mutacje? sposób określania wyjściowej struktury genetycznej?
-metody analizy struktury genetycznej
-wiek???

ZROBIONE:
-płeć osobników
  -plec jako cecha osobnikow (w pop poczatkowej przydzielana losowo)
  -dodaje Nm, Liczm 
  -samce0, samce1 -> tablica nr samcow-
-geny
	-na razie wartości początkowe losowe - random(10)
	-liczba genów jako stała ngen
	-pdb mutacji równe pmut - zwiększenie lub zmniejszenie długości allelu o 1 powtórzenie.

do zrobienia:
-sposób generowania początkowego rozkładu genów w populacjach (losowanie z "populacji zewnętrznej"?)
-statystyki:
 -średnia liczba alleli (allelic richness)
 -
}


Uses Crt;
const
 ngen=5; //liczba genów
 k=6; //liczba powierzchni
 aL=0;
 bL=-10;
 ar=-0.005;
 br=1;
 as=0.005;
 bs=-1;
 ae=0.01;
 be=-2;
 cc=0;
 czas=500;
 lpowt=1;
 pmut=0.001;
 skos=0.5; //pdb wydłużenia motywu
 maxNAllel=1000;

Type stanosobnika=record
          NRsubpop:longint;
          NR:longint;
          plec: integer; //0=samica, 1=samiec
          geny: array[1..ngen,0..1] of longint; //pierwsza cyfra -> nr. genu; 0..1 -> 0 od matki, 1 od ojca; wartość-> liczba powtórzeń motywu 
          end;
 stanpopulacji=array[1..100000] of stanosobnika;
 nrSamcow=array[1..100000] of longint;
 liczebnosc=array[1..k] of longint;
 tablica=array[1..k, 1..k] of real;
 GenDict=object
          labels: array[1..maxNAllel] of longint; //wartości alleli
          values: array[1..maxNAllel] of longint; //liczba wystąpień allelu
          max: longint;
          constructor init; 
          procedure add(lab:longint);
          procedure clear;
          function getValue(lab:longint):longint;
          end;
 GenDictArray=array[1..k,1..ngen] of GenDict;
 constructor GenDict.init;
          begin 
          max:=0;
          end;
 procedure GenDict.add(lab:longint);
          var stop: boolean;
              i:longint;
          begin 
           stop:=false;
           for i:=1 to max do
           begin
            if labels[i]=lab then
             begin
             values[i]:=values[i]+1;
             stop:=true;
             break;
             end;
           end;
           if stop=false then
           begin
            labels[max+1]:=lab;
            values[max+1]:=1;
            max:=max+1;
           end;
          end;
 procedure GenDict.clear;
          begin
           max:=0;
          end;
  function GenDict.getValue(lab:longint):longint;
          var val, i:longint;
          begin
           val:=0;
           for i:=1 to max do
           begin 
            if labels[i]=lab then val:=values[i];
           end;
          getValue:=val;
          end;

Var i, j, l, os, pot, t, powt, ll, dc, lm, s: longint;
    N0, N, Nm, Licz, Liczm, ost, lmigr : liczebnosc;
    pr, ps, pe, sum  : real; //zmienna suma nie jest uzywana
    polepow : array[1..k] of real;
    minodl,tabc : tablica;
    osob,potom,ojciec : stanosobnika;
    POP0,POP1,POPmigr : array[1..k] of stanpopulacji;
    samce0, samce1: array[1..k] of nrSamcow; 
    znak : char;
    INFO, INFO1,INFO2 : text;
    wiersz : string[k];
    allels : GenDictArray;
    avgAllelRichness: array[1..k] of real;

function normal(mi,sigma:real):real; //losuje liczbe z rozkladu normalnego
Var alfa,r1,r2:real;
begin
  alfa:=2*pi*random;
  r1:=random;
  r2:=random;
  if r2<r1 then r2:=1-r2;
  normal:=sigma*cos(alfa)*sqrt(-2*ln(sqr(r2)))+mi;
end;

Function Lpotom(n:longint;w:real;a,b:real):longint; //losuje liczbe potomkow z przesunietego o 1 rozkladu poissona o sredniej lambda
var L:longint;
    x,y,z,lambda,LOT:real; //"z" nie jest uzywany 
  begin
 if a*N/w+b>50 then lambda:=a*N/w+b else lambda:=ln(1+exp(a*N/w+b));
  if lambda<50 then
    begin
    L:=0;
    x:=exp(-lambda);
    y:=exp(-lambda);
    LOT:=random;
    while LOT>x do
      begin
      L:=L+1;
      y:=y*lambda/L;
      x:=x+y;
      end;
    end
  else L:=round(normal(lambda,sqrt(lambda)));
  LPOTOM:=L+1;
  end;

function POPDOC(k,j:longint; tabc:tablica):integer; //losuje populacje docelowa 
  var i:longint;
    czsum, los : real;
  begin
  sum:=0;
  los:=random;
  czsum:=0;
  i:=1;
  while los>czsum do
    begin
    czsum:=czsum+tabc[i,j];
    i:=i+1;
    end;
  POPDOC:=i-1;
  end;

begin
{Wczytanie wartosci parametrow modelu: k, N0, aL, bL, ar, br, as, bs, ae, be, cc, czas, lpowt}
{N0[i]] - liczebnosc poczatkowa i-tej populacji,
aL[i,j], bL[i] - parametry sluzace do wyliczania sredniej jednorazowej
liczby potomkow jednego osobnika z i-tej populacji,
ar[i,j], br[i] - parametry sluzace do wyliczania prawdopodobienstwa
rozrodu osobnika z i-tej populacji,
as[i,j], bs[i] - parametry sluzace do wyliczania prawdopodobienstwa
smierci osobnika z i-tej populacji
czas, lpowt}

for i:=1 to k do N0[i]:=10;

{wczytanie danych z pliku POLEPOW.TXT}
  assign(INFO1,'POLEPOW.TXT');
  reset(INFO1);
  readln(INFO1,wiersz);
  for i:=1 to k do readln(INFO1,znak,POLEPOW[i]);
  close(INFO1);
  writeln('wczytane pola powierzchni platow srodowiska');
  for i:=1 to k do writeln(i,' ',polepow[i]:7:2,' ');

{wczytanie minimalnych odleglosci miedzy platami środowiska}
  assign(INFO2,'MINODL.TXT');
  reset(INFO2);
  readln(INFO2,wiersz);
  for i:=1 to k do
    begin
    read(INFO2,znak);
    for j:=1 to k do 
     begin
     read(INFO2,MINODL[i,j]);
     write(' ',MINODL[i,j]:7:2,' ');
     end;
    writeln();
    readln(INFO2);
    end;
  close(INFO2);

{utworzenie tablicy tabc[i,j] na podstawie cc i minimalnych odległosci}
  for j:=1 to k do
    begin
    sum:=0;
    for i:=1 to k do if j<>i then sum:=sum+1/(exp(cc*ln(MINODL[i,j]))+1);
    for i:=1 to k do
      begin
      if i=j then tabc[i,j]:=0 else tabc[i,j]:=1/(exp(cc*ln(MINODL[i,j]))+1)/sum;
      end;
    end;
{utworzenie naglowka pliku wyjsciowego}
assign(INFO,'infodyn.txt');
 rewrite(INFO);
 append(INFO);
 write(INFO,'Powtorzenie czas ');
 for i:=1 to k do write(INFO,'N',i,' ');
 for i:=1 to k do write(INFO,'Nm',i,' ');
 for i:=1 to k do write(INFO,'Zag',i,' ');
 for i:=1 to k do write(INFO,'avgAllelRichness',i,' ');
  writeln(INFO);

{Utworzenie tablicy allels}
 for i:=1 to k do
  begin
  for j:=1 to ngen do
   begin
   allels[i,j].max:=0;
   for l:=1 to maxNAllel do
    begin
    allels[i,j].labels[l]:=0;
    allels[i,j].labels[l]:=0;    
    end;
   end;
  end;

{poczatek symulacji}
  randomize;
  for powt:=1 to lpowt do
    begin

    {utworzenie populacji początkowej}
    for i:=1 to k do 
      begin
      {wyzerowanie wskaźników}
      lmigr[i]:=0;
      avgAllelRichness[i]:=0;
      for j:=1 to ngen do allels[i,j].clear;
      Nm[i]:=0;

      for os:=1 to N0[i] do
        begin
        osob.NRsubpop:=i;
        osob.nr:=os;
        osob.plec:=random(2);
	for j:=1 to ngen do
         begin
         for l:=0 to 1 do
          begin
          osob.geny[j,l]:=random(10); //DO ZMIANY!!!
          allels[i,j].add(osob.geny[j,l])
          end;
         end;
        if osob.plec=1 then 
         begin
         Nm[i]:=Nm[i]+1;
         samce0[i][Nm[i]]:=os; //Zapisuje w tablicy nr samca w subpopulacji
         end;
        POP0[i][os]:=osob;
        end;
      end;
    for i:=1 to k do ost[i]:=N0[i];
    for i:=1 to k do N[i]:=N0[i];
    t:=0;

    for i:=1 to k do 
     begin
     avgAllelRichness[i]:=0;
     for j:=1 to ngen do
      begin
      avgAllelRichness[i]:=avgAllelRichness[i]+allels[i,j].max;
      end;
     avgAllelRichness[i]:=avgAllelRichness[i]/ngen;
     end;  

    write(INFO,powt,' ',t,' ');
    for i:=1 to k do write(INFO,N0[i],' ');
    for i:=1 to k do write(INFO,Nm[i],' ');
    for i:=1 to k do write(INFO,N0[i]/polepow[i]:7:5,' ');
    for i:=1 to k do write(INFO,avgAllelRichness[i]:7:5,' ');
    writeln(INFO);
    ll:=0; //parametr przerywajacy symulacje, jesli ktoras z populacji jest zbyt liczna
    for i:=1 to k do if N0[i]>100000 then ll:=1;

    {kolejne kroki czasowe}
    while ((t<czas)and(ll=0)) do
      begin
      t:=t+1;

      {wyzerowanie liczników}
      for i:=1 to k do
       begin 
       lmigr[i]:=0;
       avgAllelRichness[i]:=0;
       for j:=1 to ngen do allels[i,j].clear;
       end;

      {przeglądanie populacji}
      for i:=1 to k do
        begin
        Licz[i]:=0;
        Liczm[i]:=0;
        pr:=1/(1+exp(-(ar*N[i]/polepow[i]+br)));
        ps:=1/(1+exp(-(as*N[i]/polepow[i]+bs)));
        pe:=1/(1+exp(-(ae*N[i]/POLEPOW[i]+be)));
        os:=0;
        while os<N[i] do
          begin
          os:=os+1;
          osob:=POP0[i][os];
          if ((osob.plec=0) and (random<pr) and (Nm[i]>0)) then for pot:=1 to Lpotom(n[i],polepow[i],aL,bL) do
            begin
            ojciec:=POP0[i][samce0[i][(random(Nm[i])+1)]];
            {if ojciec.plec=0 then //czy wybieranie samca dobrze działa?
             begin
             writeln('ojciec to samica',' t=',t);
             Exit;
             end;}
            potom.NRsubpop:=i;
            potom.nr:=ost[i]+1;
            potom.plec:=random(2);
            for j:=1 to ngen do
             begin
              potom.geny[j,0]:=osob.geny[j,(random(2))];
              potom.geny[j,1]:=ojciec.geny[j,(random(2))];
              for l:=0 to 1 do if random<pmut then
               begin
               if random<skos then potom.geny[j,l]:=potom.geny[j,l]+1
               else potom.geny[j,l]:=potom.geny[j,l]-1;
               end;
             end; 
            ost[i]:=ost[i]+1;
            Licz[i]:=Licz[i]+1;
            POP1[i][Licz[i]]:=potom;
            for j:=1 to ngen do
             begin
             allels[i,j].add(potom.geny[j,0]);
             allels[i,j].add(potom.geny[j,1]);
             end;
            if potom.plec=1 then
             begin
             Liczm[i]:=Liczm[i]+1;
             samce1[i][Liczm[i]]:=Licz[i];   
             end; 
            end;
          if Licz[i]=100000 then os:=N[i]+1 else
            begin
            if random<ps then else
              begin
              if random<pe then
                begin
                dc:=POPDOC(k, i, tabc);
                POPmigr[dc][lmigr[dc]+1]:=osob;
                lmigr[dc]:=lmigr[dc]+1;
                end
              else
                begin
                POP1[i][Licz[i]+1]:=osob;
                Licz[i]:=Licz[i]+1;
                for j:=1 to ngen do
                 begin
                 allels[i,j].add(osob.geny[j,0]);
                 allels[i,j].add(osob.geny[j,1]);
                 end;
                if osob.plec=1 then
                 begin
                 Liczm[i]:=Liczm[i]+1;
                 samce1[i][Liczm[i]]:=Licz[i];   
                 end;
                end;
              end;
            end;
          if Licz[i]=100000 then ll:=1;
          end;
        if Licz[i]=100000 then ll:=1;
        N[i]:=Licz[i];
        Nm[i]:=Liczm[i];
        end;
      for i:=1 to k do for lm:=1 to lmigr[i] do
        begin
        N[i]:=N[i]+1;
        POP1[i][N[i]]:=POPmigr[i][lm]; 
        osob:=POPmigr[i][lm];
        for j:=1 to ngen do
         begin
         allels[i,j].add(osob.geny[j,0]);
         allels[i,j].add(osob.geny[j,1]);
         end;
        if osob.plec=1 then
         begin
         Nm[i]:=Nm[i]+1;
         samce1[i][Nm[i]]:=N[i];   
         end;
        end;

{statystyki}
      for i:=1 to k do 
       begin
       avgAllelRichness[i]:=0;
       for j:=1 to ngen do avgAllelRichness[i]:=avgAllelRichness[i]+allels[i,j].max;
       avgAllelRichness[i]:=avgAllelRichness[i]/ngen;
       end;   

      write(INFO,powt,' ',t,' ');
      for i:=1 to k do write(INFO,N[i],' ');
      for i:=1 to k do write(INFO,Nm[i],' ');
      for i:=1 to k do write(INFO,N[i]/polepow[i]:7:5,' ');
      for i:=1 to k do write(INFO,avgAllelRichness[i]:7:5,' ');
      writeln(INFO);
      for i:=1 to k do for os:=1 to N[i] do POP0[i][os]:=POP1[i][os];
      for i:=1 to k do for s:=1 to Nm[i] do samce0[i][s]:=samce1[i][s];
      end;
    end;
  close(INFO);
  writeln('Koniec');
  repeat until Keypressed;
  znak:=readkey;
end.

