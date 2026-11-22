# DRUtES GUI — návod pro studenty

Webová aplikace: **[https://drutes.org/~drutesgui/](https://drutes.org/~drutesgui/)**

DRUtES GUI umožňuje připravit jednorozměrnou simulaci proudění vody v půdě
(Richardsova rovnice) nebo přenosu tepla, spustit výpočet a prohlédnout výsledky.
Není potřeba instalovat program na vlastní počítač. Výpočet běží na serveru;
webový prohlížeč slouží k nastavení a zobrazení výsledků.

Názvy tlačítek a sekcí jsou v aplikaci anglicky. Hodnoty fyzikálních parametrů,
počáteční a okrajové podmínky přebírejte ze zadání cvičení.

## 1. Přihlášení a projekt

1. Otevřete adresu aplikace a klikněte na **Continue with Google**.
2. Vyberte svůj Google účet.
3. Do **Project name** zadejte název projektu, například `cviceni-01`, a klikněte
   na **Create project**. Název začněte písmenem nebo číslicí; nejjednodušší je
   používat písmena bez diakritiky, číslice a pomlčky.
4. U nového projektu otevřete **Reconfigure project**.

Projekt obsahuje vlastní kopii výchozí konfigurace. Při dalším přihlášení stejným
účtem jej vyberte v **Existing projects** a klikněte na **Open project**.
**Switch project** slouží k přepnutí mezi projekty, **Log out** k odhlášení.

## 2. Základní nastavení — `global.conf`

Na stránce **Global configuration** nastavte:

- **Problem & solver:** Richardsovu rovnici nebo přenos tepla podle zadání.
- **Time settings:** dobu simulace a nastavení časového kroku.
- **Jednotky:** pro Richardsovu rovnici čas `sec`, `min`, `hrs` nebo `day`
  a délku `mm`, `cm` nebo `m`. Pro přenos tepla jsou jednotky pevně sekundy a metry.
- **Observation times:** časy, ve kterých chcete uložit prostorový profil řešení.
  Musí ležet uvnitř doby simulace.
- **Observation points:** souřadnice míst, ve kterých chcete sledovat vývoj
  veličin v čase.

Klikněte na tlačítko pro uložení globální konfigurace. Pod formulářem se objeví
nabídka dalších konfiguračních stránek. Pokud aplikace hlásí chybu, opravte ji
před pokračováním. Všechny hodnoty zadávejte v příslušných zvolených jednotkách;
po změně jednotek znovu zkontrolujte i ostatní parametry modelu.

## 3. Geometrie, vrstvy a síť — `drumesh1d.conf`

Otevřete editor sítě. Na stránce **1D mesh configuration** nastavte:

1. Rozsah výpočetní oblasti v **Geometry**.
2. Intervaly a krok sítě `dx` v **Mesh density intervals**. Menší `dx` znamená
   jemnější síť, více elementů a obvykle delší výpočet.
3. Počet, polohu a rozhraní půdních vrstev v **Soil layers**.

Vedle půdního profilu se zobrazuje náčrt sítě a polohy pozorovacích bodů.
Náčrt se při úpravách aktualizuje. Zkontrolujte, že vrstvy pokrývají celou
oblast, síť odpovídá její délce a pozorovací body neleží mimo profil.

Červené zprávy nahoře vysvětlují nesoulad. Modrá zpráva potvrzuje vzájemné
zarovnání pozorovacích bodů, vrstev a sítě. Poté klikněte na
**Save mesh configuration** a pokračujte přes nabídku pod formulářem.

## 4. Numerický řešič — `solver.conf`

Na stránce **Solver configuration** nastavte řešič a jeho parametry podle
zadání cvičení. Pokud si nejste jistí významem numerických parametrů, ponechte
hodnoty doporučené vyučujícím.

Volitelně lze zapnout záznam času řešení lineárních soustav. Po výpočtu pak bude
k dispozici zobrazení souboru `out/solver.time`.

Klikněte na **Save solver configuration**.

## 5. Fyzikální parametry modelu

### Proudění vody — `matrix.conf`

Počet vrstev se přebírá z editoru sítě. Zde jej přímo neměníte.
Pro každou vrstvu vyplňte hydraulické parametry podle zadání.

V **Initial condition setup** zvolte typ počáteční podmínky a její hodnotu:

- konstantní tlaková výška (`hpres`),
- konstantní celková hydraulická výška (`H_tot`),
- konstantní objemová vlhkost (`theta`).

V **Boundary conditions setup** nastavte zvlášť **bottom boundary** (dolní
hranici) a **top boundary** (horní hranici). Dostupné typy zahrnují Dirichletovu
a Neumannovu podmínku, zvláštní výškově závislou Dirichletovu podmínku, volné
odvodnění, výronovou plochu a atmosférickou podmínku. Typ, hodnotu a znaménko toku
volte podle fyzikálního významu a zadání; jednotlivé typy nejsou zaměnitelné.

- Pro časově nezávislou podmínku zvolte časovou závislost **No** a zadejte hodnotu.
- Pro časově proměnnou podmínku zvolte **Yes** a nahrajte textový datový soubor.
- Atmosférická podmínka vyžaduje časově závislá data.

U běžných časově závislých podmínek očekává aplikace dvě číselné kolony
(čas a hodnota); u atmosférické podmínky tři. Přesný význam veličin a jejich
jednotky musí odpovídat zadání příslušné okrajové podmínky. Při nesprávné
struktuře se zobrazí **incorrect data file uploaded**. Opravte soubor a nahrajte
jej znovu. Aplikace ukládá data pod názvem `101.bc` pro dolní a `102.bc` pro horní
hranici.

Při zapnutí odběru vody kořeny se zpřístupní také editor `root4uptake.conf`.
Vyplňte jeho parametry podle zadání. Richardsovu konfiguraci uložte tlačítkem
**Save Richards configuration**.

### Přenos tepla — `heat.conf`

Vyplňte tepelné parametry jednotlivých vrstev a v **Initial temperature setup**
jejich počáteční teploty. Na horní a dolní hranici jsou k dispozici Dirichletova
a Neumannova podmínka; časovou závislost a nahrávání dat nastavíte obdobně jako
u proudění vody.

Pokud zapnete výpočet konvekce pomocí Richardsovy rovnice, pole pro ruční
zadání konvekce zmizí a musíte nakonfigurovat také `matrix.conf`.
Výpočet pak řeší vodu i teplo. Po vypnutí propojení se ruční pole opět objeví.

Klikněte na **Save heat configuration**.

## 6. Uložení celé konfigurace a spuštění

Jednotlivé stránky průběžně ukládejte jejich tlačítkem **Save … configuration**.
Po uložení se pod formulářem nabízí přechod na další stránku nebo návrat na již
navštívenou stránku. Označení **Return …** znamená možnost vrátit se do editoru;
samo o sobě nepotvrzuje správnost všech jeho parametrů.

Po návštěvě požadovaných editorů se zpřístupní **Save all configurations**.
Před jeho použitím ověřte, že jste uložili změny na každé jednotlivé stránce.
Toto tlačítko dokončí konfiguraci a otevře stránku pro spuštění modelu.

1. Případné opravy proveďte přes **Re-edit …**.
2. Klikněte na **Run model**.
3. Sledujte výpis v **Model terminal**. Během výpočtu se obnovuje přibližně
   každých pět sekund a ukazuje omezený náhled posledních řádků.
4. Chcete-li běžící výpočet zastavit, použijte **Kill simulation**.
5. Počkejte na **DRUtES finished successfully.** Při chybě si přečtěte konec
   výpisu a log výpočtu.

Opětovné spuštění stejného projektu může přepsat jeho předchozí výsledky.
Před dalším během si potřebná data stáhněte. Pro porovnávání variant je vhodné
použít samostatné projekty s výstižnými názvy.

## 7. Zobrazení výsledků

### Časové řady v pozorovacích bodech

1. V **Observation points over time** vyberte jeden nebo více bodů.
2. Pro každý vybraný bod se objeví vlastní nabídka veličin. Vyberte, co chcete
   vykreslit; jeden graf může obsahovat několik bodů i veličin.
3. Nastavte **From** a **To**, tedy zobrazovaný časový interval.
4. Nastavte **Maximum points in the entire chart** a klikněte na **Plot**.

Pro Richardsovu rovnici můžete podle dostupného výstupu sledovat tlakovou
výšku, vlhkost, stupeň nasycení, Darcyho tok a kumulativní tok. Pro přenos tepla
se zobrazují příslušné teplotní a tokové veličiny. Popisky os obsahují jednotky.
Při společném zobrazení různých veličin věnujte pozornost jejich odlišným
jednotkám a velikostem.

### Prostorové profily

V **Profiles at observation times** vyberte veličinu či více veličin a uložené
časy. Profily zobrazují hodnotu podél souřadnice výpočetní oblasti. K dispozici
může být počáteční stav, zadané pozorovací časy i konečný stav. Nastavte rozsah
souřadnice a klikněte na **Plot**.

Při propojeném výpočtu vody a tepla jsou dostupné výsledky obou rovnic.

### Velké grafy a přiblížení

Graf má limit nejvýše **10 000 bodů celkem pro všechny křivky**. U delších řad
aplikace ponechává krajní body a minima a maxima v jednotlivých úsecích.
Zobrazený graf je proto náhled, nikoli vždy úplná datová řada.

Přiblížení grafu pracuje pouze s tímto náhledem. Pro podrobnější pohled zúžte
**From/To** a znovu klikněte na **Plot**. Na telefonu nebo pomalejším zařízení
snižte limit bodů nebo počet současně vykreslených křivek. Stažené datové
soubory obsahují kompletní původní data.

## 8. Log a stahování souborů

- Tlačítkem pro zobrazení `out/DRUtES.log` otevřete přehled průběhu výpočtu.
- U výsledků lze stáhnout jednotlivé datové soubory.
- **Download all outputs as ZIP** stáhne všechny výstupy projektu.
- V postranním panelu nejprve klikněte na **Prepare configuration ZIP** a poté
  na **Download project configuration**, chcete-li uložit konfiguraci projektu.

Pro odevzdání cvičení postupujte podle požadavků vyučujícího. Obvykle se hodí
uložit konfiguraci, výstupní data a krátké vyhodnocení grafů s uvedenými jednotkami.

## 9. Když se něco nedaří

| Problém | Co zkontrolovat |
| --- | --- |
| Nelze pokračovat z globální konfigurace | Pozorovací časy musí ležet v době simulace. Přečtěte chybovou zprávu. |
| Červená zpráva v editoru sítě | Rozsah domény, rozhraní vrstev, intervaly sítě a souřadnice pozorovacích bodů. |
| Chyba nahraného souboru | Počet číselných kolon a formát dat pro zvolený typ okrajové podmínky. |
| Není dostupné Save all configurations | Projděte všechny požadované editory, včetně případné vody při propojení s teplem nebo odběru kořeny. Jednotlivé změny uložte. |
| Graf se nezobrazil | Vyberte bod, jeho veličiny a interval obsahující data; klikněte na Plot. |
| Prohlížeč při grafu padá | Snižte limit bodů, počet křivek a rozsah; zkuste počítač místo telefonu. |
| Výpočet skončil chybou | Stáhněte log, poznamenejte název projektu a předejte zprávu vyučujícímu. |
| Přihlášení nebo stránka nefunguje | Obnovte stránku a ověřte připojení. Při opakované chybě kontaktujte vyučujícího. |

Po skončení práce stáhněte potřebné soubory a na sdíleném počítači se odhlaste.
