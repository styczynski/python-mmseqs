//  main.cpp
//
//  Created by Martin Steinegger on 26.11.12.
//  Copyright (c) 2012 -. All rights reserved.
//

#include <fstream>
#include <iostream>
#include "msaFilter.h"
#include "multipleAlignment.h"
#include "pSSMCalculator.h"
#include "parameters.h"
#include "profileStates.h"
#include "sequence.h"
#include "stripedSmithWaterman.h"
#include "substitutionMatrix.h"

const char* binary_name = "test_profilestates";

int main(int, const char**) {
  Parameters& par = Parameters::getInstance();
  SubstitutionMatrix subMat(par.scoringMatrixFile.aminoacids, 2.0, 0.0);

  std::cout << "Subustitution matrix:";
  SubstitutionMatrix::print(subMat.subMatrix, subMat.num2aa,
                            subMat.alphabetSize);
  //   BaseMatrix::print(subMat.subMatrix, subMat.alphabetSize);
  const char* seqs[1001];
  int counter = 0;

  seqs[counter++] =
      "QDELTAGPCATVHVITVQMAKSGELQAIAPEVAQSLAEFFAVLADPNRLRLLSLLARSELCVGDLAQAIGVS"
      "ESAVSHQLRSLRNLRLVSYRKQGRHVYYQLQDHHIVALYQNALDHLQECR";
  seqs[counter++] =
      "-----------------"
      "EIAPLAELQAIAPEVAQSLAEFFAVLADPNRLRLLSLLARSELCVGDLAQAIGVSESAVSHQLRSLRNLRLV"
      "SYRKQGRHVYYQLQDHHIVALYQNALDHLQECR";
  seqs[counter++] =
      "---------------------"
      "AAELQAIAPEVAQSLAEFFAVLADPNRLRLLSLLARSELCVGDLAQAIGVSESAVSHQLRSLRNLRLVSYRK"
      "-----------------------------";
  seqs[counter++] =
      "------------------------"
      "LQGLELEKAQKMAEFFSLLGDANRLRILSLLAKQELCVCDLADDLGMSESAVSHQLRTLRALRLVKYQKQGR"
      "RVFYRLADHHVLDLYYAVSEHLEE--";
  seqs[counter++] =
      "---------------------"
      "SQEIQVLSSEKAQRMAEFFSFLGDANRLRILSLLAEKEFCVSDLAARLDMSESAVSHQLRNLRAMRLVNYRK"
      "QGRRVFYRLHDNHVLHLYQAVAEHLDE--";
  seqs[counter++] =
      "-----------------------"
      "QVRQVQPEVAQQMAEFFSALADPSRLRLMSALARQELCVCDLAAAMKVSESAVSHQLRILRSQRLVKYRRVG"
      "RNVYYSLADNHVMNLYREVADHLQE--";
  seqs[counter++] =
      "--------------------------------"
      "AQQMAEFFSALADPSRLRLMSALARQELCVCDLAAAMKVSESAVSHQLRILRSQRLVKYRRVGRNVYYSLAD"
      "NHVMNLYREVADHLQE--";
  seqs[counter++] =
      "-------------------------"
      "EVLSTEKSQRMAEFFSFLGDANRLRILSFLATKELCVSDIATLLEMSESAVSHQLRNLRAMRLVSYRKQGRH"
      "VFYRLHDNHILELYQAVAEHLDE--";
  seqs[counter++] =
      "------------------------------"
      "EKAQRMAQFFGLLADTNRLRIVDLLAQGEFCVRDIAVALEMSESAVSHQLRMLKALRLVRFRRQGRHIFYQL"
      "LDHHVLTLYKAAAEHLDE--";
  seqs[counter++] =
      "-----------------------"
      "QTQVLNSQKAQRMAEFFSLLGDANRLRLLSVLAAQELCVCDLAATLEMSESAVSHQLRALRALRLVSYRKQG"
      "RQVFYSLLDRHVLELYRAVAEHLDE--";
  seqs[counter++] =
      "--------------------------------"
      "AQQMAEFFAVLADPNRLRLISALASQELCVCDLAALMKMTESAVSHQLRLLKAMRLVSYRREGRNIYYSLAD"
      "NHVISLYREVAVHLDE--";
  seqs[counter++] =
      "---------------------------------------"
      "FAALGDPTRFRIIAALQVQELCVGDLAAAIGLSQSAVSHQLRALRDLGLVRSRREGRLVYYALDDEHVVTLV"
      "AQALDHVREER";
  seqs[counter++] =
      "-------------------------------"
      "VTRQMAEFFKSLSDPTRLRIVQALLEEELCVCDISAIVDISISAISHQLRLLRSMHIVKFRKQGKMVYYSLE"
      "DEHISRMLEIALEHLNE--";
  seqs[counter++] =
      "------------------------------"
      "DITNRLAETFKVLGDPTRLKILLAVSLDELCVCDIASLLGTTKSAVSHQLRLLRSLRVVKYRKDGRIVYYSL"
      "DDSHVGNLLSEGLDHI----";
  seqs[counter++] =
      "-----------------------------------"
      "MAETFKILADPTRVKILHALAHKELCVCDIAVTLDMKVSAVSHQLRLLKSARLVKQRREGKNVYYQLDDHHV"
      "EQLFEKTLEHIKQ--";
  seqs[counter++] =
      "--------------------------------"
      "AGRLADLFKALADPTRVRIIAALLHTELCVDDLANLLDMSQSAISHQLRLLRNLHLVQFRRSGKHAFYRLVD"
      "DHVRDLFQRSREHL----";
  seqs[counter++] =
      "------------------------------"
      "ELLYELAELFKIFGDSSRIRILSLLQKERLCVSEISTLLNLSQSAISHQLRILRQARLVRYKKIGKEVFYEL"
      "DDDHIEKIFEQGLEHIQE--";
  seqs[counter++] =
      "--------------------------------"
      "ATRLAAAFQALSDPTRVRLISALLEQELCVHDLAAVLGMSQSATSHQLRVLRALGLVRTRKEGRIVYYALDD"
      "EHIRELFQRGLEHI----";
  seqs[counter++] =
      "---------------------------------"
      "EEMASFFRMMGDPTRIRILSLLFDEELCVHTLAERLEMTHSAVSHQLALLKHARLVRSRREGRHVYYRLADE"
      "HVQKVYELAREHLEE--";
  seqs[counter++] =
      "-DEFKKGSDLNCHVMNIQDCKEVELKSVTEETIFDLSEFFKVFSDSTRIKILSSLLVSEMCVCDLAAVLGTS"
      "QSAISHQLRLLKVFRLVKSRKAGKVVYYSLSDDHVKSIIELGLAHLSE--";
  seqs[counter++] =
      "------------------------------"
      "ELIFNLADFFKTFGDSTRIKIICALMETELCVCDLANVINTSQSAVSHQLRVLRQSRLVKYRKDGKTVYYSL"
      "DDDHIKLLISQGLDHL----";
  seqs[counter++] =
      "-----------------------------------"
      "MAEVFKALNDPTRLKIINILIVSELCVNDIANLLEISQPAISHHLKELRQLKLIKYHKKGRSVFYSLDDEHI"
      "HPLFQQCLEHVNE--";
  seqs[counter++] =
      "-----------------------------------"
      "LAELFKILGDPTRLKIVELLLENEMCVNHIAETMGMGQSAISHQLRVLRQARLVTYRKDGKTAYYSLNDNHV"
      "ECLVRMGMEHV----";
  seqs[counter++] =
      "---------------------------------"
      "EDLALLFKMFADPTRLKVLKALFEREMCVGDLAVLLKMTHSAVSHQLASLKKTRLVRSRKDGKVVYYSLDDD"
      "HIEEIFQKALDHVRE--";
  seqs[counter++] =
      "-------------------------------------------"
      "ADQTRLRILCLLRDREVCVHDIVEALDMSQSAISHQLRVLRDARLVSHRREGRHVYYRLADDHVREMLENAL"
      "SH-----";
  seqs[counter++] =
      "-----------------------------------"
      "LSMLFKMFADPTRLRIFTILSHQTVCVDDLAEILGMTQSAVSHQLASLRKMNLVRSSKVGKNAYYQLADSHV"
      "MQIFSQALDHVKE--";
  seqs[counter++] =
      "-------PCAEVTALHLRL--"
      "SPATRAVDKGALQRAAAIFRALGDPARLHLLALLAAGEQCVSQLATETGDSLPAISQRLKLLRSERLVSQRR"
      "DGKHIYYQLADQHVVQLIEAGIDHAVESR";
  seqs[counter++] =
      "-----------------------------------"
      "LSEVYKSLGDGTRLKILLALKEKESCVCDLAAALSMSQSAISHQLRVLRNVRLVKYRREGKMVYYSLDDEHI"
      "LKILQEGLNHI----";
  seqs[counter++] =
      "------------------------------"
      "EEVQKLSAIYKALGDPTRFKILFCLKQEEMCVCDISAILDMSQSAISHQLRVLRNLRIVKYRKEGKMVFYSL"
      "DDKHIFRILDEGINHIRE--";
  seqs[counter++] =
      "------------------------------"
      "ETMDAIAELFKGFADSTRVHILALLSRQELCVTDIAETVDVSQSAISHQLRILKQMHLIKFRREGKNILYSL"
      "ADDHVKTILQMGLEHVLE--";
  seqs[counter++] =
      "-------------------------------------"
      "ELLKAVGDPTRMRILCALADRELCVCDLQAVLGLSQSAVSHQLRTLRNARLVTYRREGKMAYYTLADDHVRR"
      "LLDLSLEHV----";
  seqs[counter++] =
      "-------------------------"
      "EAKADQIASPLANFFKTLSDPTRLRIILAIGTTSLSVNEISTIINMSQSSVSHQLRILRDNHLVISQRFGQH"
      "IHYQLTDQHVLTILENSLDHISE--";
  seqs[counter++] =
      "-----------------------------------"
      "LAATFRLLGDRTRVRILEALAGDELCVCDLAAVVGHSQSAVSHQLRLLRAAKLVRVRRDGKNAFYSLDDDHV"
      "RHLFRQALDHVQE--";
  seqs[counter++] =
      "-----------------------------------"
      "MAQLFKILGDPTRVRILQALSISEMCVCDIAALLEMTQSAISHQLRLLKQGRLVKYRRDGKVVYYSLNDNHV"
      "RLIFDQALSHITE--";
  seqs[counter++] =
      "--------------------------"
      "ALTEKTAKELSELFKVVSDPTRIKILWAIGGGEVCVCCISELLGMSVSAVSHQLKTLRQAHLVKARREGRNI"
      "YYSLDDHHVKILLDVLLEHMEE--";
  seqs[counter++] =
      "---------------------------------"
      "ENLGEFFKVLTDASRLKILYALGAGELCVFDLSVTIGASVSSVSHHLAALKRVRLVKGRRDGRIIYYSLDDD"
      "HVKSIIRYAREHLEE--";
  seqs[counter++] =
      "---------------------------------"
      "QKLSNMFKLFSDETRLKIICSLLKEELCVCDLCELLGLNQSQVSHQLQLLRNSKLVKFRREGKQIFYSLDDE"
      "HVELIIQMALDHILE--";
  seqs[counter++] =
      "-----------------------------------"
      "LSETFGALADSNRAKILHSLLNQELCVCDIACVVGISESAISQHLRILRTLRLVKQRKQGRMMYYSLNDNHI"
      "RQLLEICLEH-----";
  seqs[counter++] =
      "-----------------------------------"
      "LAQLFAAFGDPTRLRILTALRSGDLCVCDLTAVLGMTASAVSHQLRLLRNLRLVRSRKVGRVVYYHLDDEHV"
      "LNL------------";
  seqs[counter++] =
      "------------------------------"
      "EILGDLSDFFKVIGDGTRIRILWALDVSEMCVCDIANVLNMTKSAVSHQLRALRDADLVKFRKSGKEVLYSL"
      "SDNHVKEIFEQGLIHIQE--";
  seqs[counter++] =
      "-----------------------------------"
      "LADLFKMFADSTRLKILCILCESEMCVNDIANLISMSQSAVSHQLRILKQSKLIRGRREGKIVFYSLADSHI"
      "NTIINNGLEHIQE--";
  seqs[counter++] =
      "--------------------"
      "QSQEREVLAAPLAWRVADIFKALGDPTRVKIIALLDAGEMCVGEMCLTLGMSQPAISSQLRLLRTLGIVSVR"
      "REGKHAYYRLADEHVRHLFHQGLAH-----";
  seqs[counter++] =
      "------------------------------"
      "EYIQELSAFFKVFGEENRTRILYALSIREMCVNDLVTLLGMSQSSVSHQLQILRAHGQVKFRKEGRNVFYSL"
      "DDKHVVDVFQEALQHI----";
  seqs[counter++] =
      "------------------------------"
      "KVIYELSEFFKILSDQTRLKILVLLFEKEQNVSELQRQIGVTQSNISHQLRILRQANLVRYRKIGRNVYYRL"
      "YDEHVEIIIKYAMEHLKE--";
  seqs[counter++] =
      "------------------------------"
      "EVMFDLSEFFKVFSDSTRIKILSSLLVSEMCVCDLAAVLNTSQSAISHQLRLLKAFRLVKSRKVGKVVYYSL"
      "SDDHVKSIIELGLTHLSE--";
  seqs[counter++] =
      "-----------------------------------"
      "LAELFKIFGDATRIRILCALSEGEICVSDLAETLSMTQPAISHQLRILKNTRLVKARRDGKQIYYSLADAHV"
      "SSIIGTALEHVEE--";
  seqs[counter++] =
      "------------------------------"
      "ETLYDLADLFKVLGDSTRIKILCTLFEAEMCVCDIAAVLGMTQSAISHQLRVLKQARLVKYKRSGKVVYYSL"
      "DDEHVKHIFDQGLIHISENR";
  seqs[counter++] =
      "-----------------------------------"
      "LSDFYKVFGDPTRLKILFALESRELCVCDLAQILQMTKSAVSHQLKILRQTELVNFKKLGRSVFYRLSDAHI"
      "QGILDQGADHINE--";
  seqs[counter++] =
      "------------------------------"
      "ECVMDLADFFSIFSDSTRIRILWVLYGRELCVRDISDTLGISMSACSHQLKTLRNSGAVEARRDGKMIYYKL"
      "ADEHVEILLRTGLEHIQE--";
  seqs[counter++] =
      "-----------------------------------"
      "LADMYKALGDPSRLRIVMALSQGEMCVCDLAAYLEISESAVSHQLRRLRSLALVKNRRDGKILYYSLDDDHV"
      "SSLVALGLEHVRE--";
  seqs[counter++] =
      "-------------------------------------"
      "DLFKTLGDPSRLRIIEVLSQEELSVDDLANQVGLTQSAASHQLRRLKLDRVVKYRKEGKFIYYSLADQHLLY"
      "LFNIARDHVQE--";
  seqs[counter++] =
      "-----------------------------"
      "PDQIESLSNFYKIMGDPTRLMLLMALEAGELCASDLANVTNMSRSAVSHQLKTLKQACLVRSRRDGKTIFYE"
      "LDDEHIYSVLKVAFEHIQE--";
  seqs[counter++] =
      "-----------------------------------"
      "LSDFFKVLGDSTRARIISALDINEMCVCDLAVLLNMTKSAISHQLRSLKEANLVKFRKEGKVVFYSLTDDHV"
      "KDIFEKGLEHIRE--";
  seqs[counter++] =
      "-----------------------------------"
      "LSEFFKVFSDSTRVKILSALLISEMCVCDLAALLQVTQSAISHQLRLLKAFRLVKSRKEGKVVYYSLNDDHV"
      "KSILELGLLHLSEAK";
  seqs[counter++] =
      "-----------------------------"
      "PDIIDDLSELFKILGDQTRSKILFVLEQGEFCVSDISEAVGMTKYAVSHQLRTLKQAKLVKCRREGKEVIYS"
      "LDDDHVSTLFSCALAHVTE--";
  seqs[counter++] =
      "------------------------------------"
      "AELFKVFADSTRVKIINVLLENKLCVGDIAALVGGTQSAISHQLRILKSAKLVKYTKIGKTVYYELSDDHVK"
      "KLFSVGKEHINE--";
  seqs[counter++] =
      "------------------------------"
      "ETMSDLAAIFKLMGEPVRITILHALSIRDLCVCDLAELLGMSHSAVSHQLRLLRTARMVRFEKQGRKAIYSL"
      "NDRHVETIMQTALAHMQ---";
  seqs[counter++] =
      "---------------------------------"
      "QELANNFKVLGDPTRLRILQALMHGERCVRELADGIQMEQSAVSHQLRTLRDAGLVNFRRDGKVVYYSLADA"
      "HVFTLLSVGIEHVAE--";
  seqs[counter++] =
      "-----------------------------------"
      "LSELFKILGDYTRIKIIYSLSKKELCVCDISEVVQMSQSAISHQLRILKAARLVKFRREGKSVYYSLDDEHI"
      "DRLFNAGLEHVE---";
  seqs[counter++] =
      "------------------------------"
      "EIINELSEFFKVFSDTTRLRILEVLLNEETSVGVIANKINVSNSAVSHQLSYLRSTNLVKTRKEGQVIYYSI"
      "ADNHIKVIIEYGLEHIKELK";
  seqs[counter++] =
      "-----------------------------------"
      "MSNFFKAISDPTRLRILQAVRQKTICVGDLAIALQMTKSAISHQLRYLRDCQLVKGEKKGKMTYYELADDHV"
      "AAVLSLTLKHLKE--";
  seqs[counter++] =
      "------------------------------"
      "DVVASLSELFKALGDPTRVKILSCLQISDMNVGDIADKLGMTTSAVSHQLRVLRAIKLVKGTKEGKEVRYSL"
      "DDDHVTLIMQCGLTHVNE--";
  seqs[counter++] =
      "-----------------------------------"
      "LAELFKVLGDYTRIKIIYALLKKELCVCDIAELLDMSQSSISHQLRTLKAARLVKFRKEGKVVYYSLDDEHI"
      "EHILNASLEHVE---";
  seqs[counter++] =
      "----------------------------------"
      "SIADFFAVFGDRTRIKILLALDQSPMCVCDLAVLLDMTKSVISHQLSSLKKINLVSSHKEGKHSYYALADDH"
      "IKKIIEMAVEHLEE--";
  seqs[counter++] =
      "------------------------"
      "VKSIDADTAQHLADLFKTLGDPTRIKILSLLTKTEMRVYDIADSLTMGQSAISHQLRVLRSARLVKFRRDGK"
      "EVLYSIDDDHVMKLLGQGLEHVQ---";
  seqs[counter++] =
      "---------------------------------"
      "KDLADTFKLLSDFTRVKILYVLSLSELCVQDISELTGVSQSAVSHQLRILRNSRLVSWKKTGKQVFYSLNSD"
      "AVRALIEKGMEHV----";
  seqs[counter++] =
      "-----------------------------------"
      "LAELFKVFGDTTRVKILFALFTAEMCVCDLTALLGLTQSAVSHQLRVLKQARLVKYRKDGKMVYYSLDDDHV"
      "KQIFDQGLAHI----";
  seqs[counter++] =
      "--------------------KSNE-"
      "QMLAPDDVDVTAEFFKALADPTRINIVNALQIHELCVTDLAEILGMTKSAVSHQLCYLRLNNLVMVKREGQR"
      "VYYALCDEHVEKVFEMAISHIKE--";
  seqs[counter++] =
      "------------------------------------"
      "ADIFRALGDPSRLRMLSLLIHDELCVTEIAEALGDNLSAVSQRLKLLKSERIVGARREGKHIFYRLSDHHV-"
      "--------------";
  seqs[counter++] =
      "------------------------------"
      "DIVAKLSDFIKVLGDGTRIKIIWILEENEMCVNDLAVALNMSQSAVSHQLKTLKTANVVKSKREGKNIFYSL"
      "SDDHVKDIFLKTLEHIQE--";
  seqs[counter++] =
      "------------------------------"
      "ETMSDLSDFFKNFGDSTRIKIVSALISGELCVADIAEVLEMSVSAISHQLRILRQAKIVRARRNGKQMYYSI"
      "DDEHVAILYSLGLEHIRE--";
  seqs[counter++] =
      "------------------------------"
      "ETFNSLSDNFKVLSDPTRLKILYALMLKEICVCDLAAVLEMTDSAVSHQLRLLRNRNLVKFRKKGKMAYYSI"
      "SDHKIV--------------";
  seqs[counter++] =
      "------------------------------------"
      "AELFKMFGDPTRLKLLAALLGQEVCVCDLSDLLGISQSAVSHQLRLLRASHLVKNRREGKSVFYSLDDEHVA"
      "TILAQGMEHV----";
  seqs[counter++] =
      "-----------------------------------"
      "LSETFKIFGNPTRLKILSLLSVEDLCVCDICEILNMSQSAVSHQLRTLRSKNLVKYTKEGKQARYSLADKHV"
      "VQILKIGIEHVLE--";
  seqs[counter++] =
      "-------------------------------------"
      "DIFKILGDPSRMRIVAALRIKELCVGDISALMEISLSGVSHQLRLLKKSRIIKSRREGKLIYYSLDDDHIES"
      "LIDIALDHVRD--";
  seqs[counter++] =
      "---------------------------"
      "LSPQIVEEASRILKAISDPTRMKILYLLFQEECSVGHMVEVLGVSQSAISHQLTHLRHLRLVKYRREGNTYF"
      "YTYDDEHVVGILHQVIQHVE---";
  seqs[counter++] =
      "-----------------------------------"
      "LAELYKVFGDSTRIRILYALLESEMCVGDMAQLLGLTPTACSHQLRVLKNSKLVRFRREGKIMFYSLADDHV"
      "RSILALGMEHILE--";
  seqs[counter++] =
      "------------------------LRPVAEEAAK-"
      "LAPIFKALSDETRVKIIYALAQAELCVCDIAELTGCTLPAVSHHLRILRNIGLAKSRKEGKFIYYNIDDHHV"
      "SQIINAAFAHLRE--";
  seqs[counter++] =
      "----------------------------------"
      "ALSELFKILADPSRLRILHALQSPERCVCDLAALLDMSQSALSHQLAILRRARLVRPVKIGKIVYYQLDDHH"
      "VDALIALAMEHVSE--";
  seqs[counter++] =
      "----------------------------"
      "APEGTRRVAEVFGVLSDPTRARIVCALSIEELCVCDVAAVAGLSVSAASHQLKRLRDRGVVDYRKEGRLAFY"
      "RLVDDHVRSLMEEGVE------";
  seqs[counter++] =
      "-----------------------------------"
      "LAELFKVFGDPTRIRILWTLNEAEMCVCDIAVLLNMTQSAISHQLRVLKQANLVKNRKEGKAVYYSLVDDHV"
      "REIFDQGLIHINE--";
  seqs[counter++] =
      "----------------------------------"
      "AVAELFKALGDGSRMSILNALLCSELCVCDLTMILKMTQSAVSHQLRVLRGAKIVKSRKEGKNVYYSLDDPH"
      "IAMLIETGFEHVRE--";
  seqs[counter++] =
      "-----------------------------------"
      "LEELFKVMGSQTRLRIIALMEANELCVEHLADTVNISVSAISHQLKGLRQLRLVKTRKQAQSVYYSLDDEHI"
      "ALLFNTARTHLSE--";
  seqs[counter++] =
      "------------------------------"
      "QVIEELAQFFKAVADASRLKLLLYLMKQEANVNELAEATGLTLSGVSQQLKLLKLMKLVKSSKQGKYVYYSL"
      "DDHHVKHIINDSLIHLSE--";
  seqs[counter++] =
      "-----------------------------------"
      "LADLFNIFGDTTRIKILHVLSKSEMCVCDISSLINMNRSSVSHQLKTLRQAKLVKYRKEGTIVYYSLSNDHV"
      "KQVFNQGLIHLIE--";
  seqs[counter++] =
      "------------------------------"
      "EVVQTLSTYFKALADDNRLKIIHALSREELCVRDVADIIGTTMQVASHHLRVLRDIGLVKNRKEGKHVFYSL"
      "KDRKTAEFVQNVIEDLED--";
  seqs[counter++] =
      "---------------------------------"
      "ERLAGVFRTFGDGTRLRILFALLRKEMCVNDLAQNLGMTVSAVSHQLRILRQGELVRTRREGKTVYYAIADH"
      "HVSLIIRSGAEHVLE--";
  seqs[counter++] =
      "-------------------------------------"
      "KLFKIYSDFTRLRIIDLLIEKEHCVQDIADSLDASQSAISHQLKLLRDLHVVKTRKQGKQVFYSLQDNHIKE"
      "IFLIGYSHATEC-";
  seqs[counter++] =
      "---------------------------------"
      "QKLSKLFKVLSDETRIKILYTISKHEVCVNDIANVLNLSQSAISHQLKTLKDANLIKSRREKQTIYYTLVDD"
      "HVHLIYNQALSHIKE--";
  seqs[counter++] =
      "------------------------------"
      "EAIQEVSKIFKMISDPTRLSILFLLQKEELSVGAIAQSLSMEQSAISHQLKTLRTSRLVKSKRAGKNMIYSL"
      "DDLHVFSILEQVLTHIEE--";
  seqs[counter++] =
      "---------------------------------------"
      "FSMFADGTRLKIMSALFVSEMCVGDLAILLNMSTSAISHQLASLKKTKLVRMRKEGKNVFYSMDDEHIEKIF"
      "QVTYLHVKE--";
  seqs[counter++] =
      "------------------------------------"
      "AALLGVISDPTRLKIMDALRLGELCVCDLAAVLSMSVSAISHQLRLLRTARLVRGRRAGKVIFYTIHDGHVE"
      "KLIDMALDHCRQC-";
  seqs[counter++] =
      "------------"
      "HTVTLEILVQGELRESIPEERQNVADIFSLLGDPARLRILIALSVRRLRVCDLVKVVGASESSVSHALRILR"
      "AHRVVDVIRRGREAHYALADSHVRALLELAIDHV----";
  seqs[counter++] =
      "-------------------------------------"
      "DFFDALGNPTRLKILALMEAGELCTCDLSAITKLSVSAISHQLRILKDRKIVTYRKDGKNVFYRLDDEHIRE"
      "ILRTALNHLSEVR";
  seqs[counter++] =
      "------------------------------------"
      "AEFFKVLGDPTRIKIVSLFENGELCVNDIVDVVDVSRTAVSHQLRILKDKRIISFRKEGQMKFYHLDDAHVE"
      "VIVLLTITHLQ---";
  seqs[counter++] =
      "---------------------------------"
      "EQLAELHKAMGDYTRIRILWYLMQKEYCVSDLAQKMEITESAISHQLRALRIARLVQSHKAGKNVIYSLQDE"
      "HIRWILEKTYDHISE--";
  seqs[counter++] =
      "---------------------------"
      "ISKELIGSAAAFFKVLGDETRMKILCTIADSEVCVNDIAEAVDMTKSAVSHQLKLLKDDDLVKSRREGKNIF"
      "YSLDDQHVMDILDIAFVHI----";
  seqs[counter++] =
      "-----------------------------------"
      "LCDIFKVLSDPTRMRIILTLVDSEMCVCDIAGAVECSQSNVSHQLRLMRQSGIVKFRKDGKSVYYSLDDDHV"
      "KTIIVQAVNHI----";
  seqs[counter++] =
      "------------------------------"
      "QAAVQFADWFKAFSDPTRVKIISALLKRELCVHDLTVLLEMGQSAVSHQLRYLRNLRIVKRRKVGKTVYYSL"
      "DDTHI---------------";
  seqs[counter++] =
      "------------------------------"
      "EILFDLADLFKVFGDSTRIRIMNVLFSGPTSVGEIAEVLDMSQSAISHQLKSLKDNNLVSSKRSGRSMYYEL"
      "ADDHVKTIFMTGLEHIKE--";
  seqs[counter++] =
      "-------------------------"
      "QQIEKNIIDSTSEFFKILADNTRLHIINLLLDREMFVNDIANALNMTNSSISHQLKKMKDNDVVKSRKDGKE"
      "VYYSLNDDHVKSIFLTTIDHI----";
  seqs[counter++] =
      "----------------------------"
      "AEEILYDVAELFKVFGDSTRIRIICALFESEMCVYDLAACLDMTQSAISHQLRILKQANLVKFRRDGKLMYY"
      "SLDDEHVKQIFDAGYKHIE---";
  seqs[counter++] =
      "----------------------------------"
      "SITQIFKAMADPTRVQILYLLSDKEYSVGEIARSLGFNQTTVSHQLRFLKNLRLVKSRREGTSIYYTQDDKH"
      "VLELLKQTIRHVE---";
  seqs[counter++] =
      "-------------------------------------"
      "QIFKALSDPNRVKIAYYLQHQELCVCDIADLLDVSVATASHHLRQLKALHIAKSRKEGKNVYYSLSDHHIQT"
      "LVAMTLEHQKEMR";
  seqs[counter++] =
      "------------------------------"
      "ETLYDVAELFKAFADTTRIKIIAILKEETLCVGAISEILNISQSAISHQLKALKNAKIVKSKREGKWIYYSL"
      "DDEHIKRIFDMGFEHI----";
  seqs[counter++] =
      "---------------"
      "TVEQETSTEHKNLPPEKVMPLADFFKVFGDSTRLKIIGILRHTSLSVCCIADCLGMEQSAISHQLKVLRNNH"
      "IVKVEKKGKQSFYSLDDLHVELLYQMGLEHIME--";
  seqs[counter++] =
      "---------------------------"
      "VTPEEATQLADLFRLLGDPTRAQLLALLEAGELCVCDLTETVEVSDTAVSHALRLLRTAGIVASRRAGRMIY"
      "YRLADVHVRMLLDLSREHL----";
  seqs[counter++] =
      "-----------------------------"
      "PAAAGKVAETLQALASPNRLRILTRLRQSPCSVTELSAAVGMEQPAVSHQLRLLRNLGLVAGDRSGRNIVYR"
      "LYDSHVASLLDEAVYHIE---";
  seqs[counter++] =
      "---------------------------"
      "IQAETLSAAAEIFKALGDVTRLRILDVLERREMKVQEIAAALGMTQSAISHQLGTLRTLRIVKARREGRSTF"
      "YSIDDGHVKQLFDLCVEHV----";
  seqs[counter++] =
      "---------------------------IAKEVS-"
      "GLSDLFKVIADETRTKIVFLLSETELCTCDLAEILRLSLPTISHHLKQLKSYRLVKSRREGKSVFYSLEDFH"
      "VVELIKLAKEHFQE--";
  seqs[counter++] =
      "-----------------------------------"
      "LAEVFRLMADPGRIRILALLEAGEVCVHDLASVSGLSESSVSQALRLMRAQRVVAGRRAGRHVFYRLDDSHV"
      "RMLLDLAITHV----";
  seqs[counter++] =
      "------------------------------"
      "EVLYDLAELFKVFGDSTRIRILYVLFETEMCVYDLSKILNVTQSAISHQLRVLKQNKLVKFRREGKNIFYSL"
      "ADEH----------------";
  seqs[counter++] =
      "------------------------------"
      "KVLNHLSGIFKVLADPTRLRIIYTLSMGELCVTDISETLEMTQSSISHQLAILKSRDLIKVRKVGRKSYYSL"
      "DDDHVLSLFEGGYEH-----";
  seqs[counter++] =
      "-------------------------------------"
      "QLFKLLGDNTRLNIITLLTRKELCVEDLMKCTGMEQSAISHQLKKLRAHHIVKAEKVGKHVWYSLEDHHVLQ"
      "LYNLADEH-----";
  seqs[counter++] =
      "--------------------------------------"
      "FFKVLADDTRIRILYALKEQEMCAGDIAVFLDMTKSAVSHQLAVMRKMHQVRARREGKNVFYSLDDQHIVDI"
      "MEEALIHM----";
  seqs[counter++] =
      "-----------------------------------"
      "MSEFFHMFDDPTRLKIIGALIISEMCVCDIAAVTGMSQSVISHHLKILRQERVIQFRRQGKMAYYSLCDDHI"
      "GGIFYQGRIHVQEER";
  seqs[counter++] =
      "---------------------------------------"
      "FKALSDPTRIRILNLLCADEHSVNDIAEILDLGQSTVSHQLRFLKNLRLVKFRREGTTLFYSKDDDHIMNLL"
      "KQAIEH-----";
  seqs[counter++] =
      "------------------------"
      "MQVLSPEQVGDVAEVFRLLGEPNRLRIVLACLETERTVGEIGEALGLSQSLTSHHLRLLRTARILRAIRHGR"
      "HVAYAIDDDHVRDVLRNMVAHLTE--";
  seqs[counter++] =
      "------------------------------------"
      "ADVFGLLGDPRRLKLLVALLEGELCVCDLAAVTGMSESATSHALRLLRAHRVVSVRRDGRMAYYRLDDAHV-"
      "--------------";
  seqs[counter++] =
      "------------------------------"
      "ETVYNMATLFSTFSDSTRLKILLCLSRAKLCSCDISAAVNISKSATSHQMRILKMTKLVKAERKGKQIFYSL"
      "SDEHVSILLQAALEHVKE--";
  seqs[counter++] =
      "---------------------------------------"
      "FQALGDSTRLQILYALMHHTLCVRDLAILVGVSESAVSHQLRLLRDRRLVRQRRSGNIIYYSLDDEHLAVLF"
      "REA--------";
  seqs[counter++] =
      "------------------------------"
      "DVFEEMSGFYKLFSDRTRLKILWALRNGPLCVCDLCAVVGMSQPAVSQQLQKLRNGRIVKSRREGKVVYYSL"
      "DDEHIEAALDMAMEHVEE--";
  seqs[counter++] =
      "-----------------------------PEALNQAAEMFRAMGDPERLRLLTMLQGGERCVGEL---"
      "VGENDSTVSARLQSLHRARLLHRRREARHIFYRLADEHVAELLNNALEHASE--";
  seqs[counter++] =
      "------------------------------------------"
      "LATPSRLYILARLQEGPCSVGDLAEAVGMEASACSHQLRLLRNLGLVTGERQGRSIVYALYDHHVAELLDQA"
      "LFHVE---";
  seqs[counter++] =
      "------------------------------------"
      "AQTFKALSDPTRIRILHLLSQGEHAVNGIAEKLNLLQSTVSHQLRFLKNLRLVKSRREGTSIYYSPEDEHVL"
      "DVLQQMIHHTQD--";
  seqs[counter++] =
      "-----------------------------------"
      "LAVIFGLLSDPGRVRILIALLEGEMCVCDLAATTGLSESGVSHALRLLRGPRVVQVRRSGRMAYYSLADSHV"
      "RMLLDLGLTHV----";
  seqs[counter++] =
      "--------------------------------"
      "AAKLAPLFKALSDETRVKIIYALAQQELCVCDIAELTGCTLPAVSHHLRLLRTMGLARCRKEGKFIYYTIDD"
      "HHVWQIINAAFEHMKE--";
  seqs[counter++] =
      "----------------------------"
      "AEHLVSAAAALFKVIGHPTRVRILLALAAEELCVCDLAQVLDATVSATSHQLRNMRAMGLVYFRTEGKLAYY"
      "RASDPVMVSLLQQGVEH-----";
  seqs[counter++] =
      "-----------------------------------"
      "LASLFKLFGDGTRVQILHALEQSEMCVCDLAVLLGVTKSAISHQLKALRLANLVKFRKEAQIVYYSLADDHV"
      "KEIIDKGFEHL----";
  seqs[counter++] =
      "---------------------------------"
      "KQVSQLYKVLSDPTRLRILLLLKEGEHNVTAISEQLGMEQSAVSHQLKLLRDSRVVKARREGKTIFYTLDDH"
      "HVIDILNQTFEHIE---";
  seqs[counter++] =
      "-----------------------------------"
      "MADLFKIFGDSTRIRILWALHESEMCVRGISKSLDLSMSAVSHQLKALKDADLVQSRRDGKNIYYSLCDEHV"
      "EILLNTALTHLKE--";
  seqs[counter++] =
      "----------------------------"
      "AEPVTEGLSTFFKALADDTRLKVIHALSQDELCVCDVANILGSTVQAASHHLRVLRNIGLAKYRKEGKRVFY"
      "SIRDRKTAGFIQKVIEDLKK--";
  seqs[counter++] =
      "------------------------------"
      "EMIEHLSEFFSMFSNPTRLRILLLLSKKEMCVGKIAEILRMDQSAVSAQLKVLRHLNLVKAKRHGRYMRYKL"
      "NNKHV---------------";
  seqs[counter++] =
      "------------------------------------"
      "ADIFKALADETRLKIFALTQENELCVCDVANIIGTSNATASHHLRHLRNLRIAKSRKEGKLVFYSLDDPHVI"
      "TLVTMAMAHGQE--";
  seqs[counter++] =
      "-----------------------------------"
      "LAELFKVFGDSTRIRILFVLFEAEVCVCDLAEALHMTQPAISHQLKILKQAKLVRSRREGKSVFYSLADGHV"
      "RTIIAQGREHIEE--";
  seqs[counter++] =
      "-----------------------------------"
      "LADVFSVLGDPGRLRLLYALRDGEVSVGALSTLTGQSDSAVSHALQLLRAHRIVRVRREGRRAYYRLDDPHV"
      "QMLLEVARSH-----";
  seqs[counter++] =
      "---------------------------------"
      "QQIAQIMRLFGDPARLRLLVLLEVDEMCVGDLAQLAQMNESATSHALRLLRAHHVVAVRREGRMAYYRIIDT"
      "HVKASLQLTLDHV----";
  seqs[counter++] =
      "-----------------------------------"
      "MAETFRALADSTRVKIVGCLLEQELCTADLAAILNYSESAVSQHLRVLRQLRLVKQRREGKLVFYSLDDDHV"
      "RVLVLVCLNHI----";
  seqs[counter++] =
      "---------------------------------------"
      "FKALGDPTRARIIYALAVSKLSVGELATGLALTQSNVSHQLTVLKQLKLVVGTRNGRNVHYQLADKHIISIF"
      "QQVAAHAEE--";
  seqs[counter++] =
      "-----------------------------------"
      "LADVFGLLGEPNRVRLLIALLNGPMCVRDLAATIEMSESAVSHALRLLRAHRVVDVHRKGRVASYELADLHV"
      "LTLLKLGLEHV----";
  seqs[counter++] =
      "---------------------------------"
      "QKVSNLFKVIADPTRIDILYTLKDSRLSVSEIKDKLNMSQSAISHQLRVLKDVNLVKDERVGKNIFYSLSDN"
      "HVYDIFNQAIDHVRE--";
  seqs[counter++] =
      "----------------------------------"
      "SIAQIFKTLSDPTRLKILYVLSKKDLCVSDISELLSMSQSSISHQLALLRHQQLIKVNRVGRMAIYSLDDDH"
      "VLSIFNQGKTH-----";
  seqs[counter++] =
      "---------------------------------------"
      "FKVLGNQSRIRILLEIADEEKCVHEISEETDQSFSNASHHLKTLRDNRLVDYRKEGRHKYYRIKDDHVLKIL"
      "QECIDH-----";
  seqs[counter++] =
      "-----------------------------------"
      "LSELFKILGDKTRINIIWTLDNREMCVCDIANVLNMTKSSISHQLAILKNAGIVKYYKSGKEVYYTLDDEHI"
      "NKLYEIGLLHIEK--";
  seqs[counter++] =
      "---------------------------------"
      "EELANFYKIFSDPTKIKILWALDISEMCGCDLAAITGVTKSAISHQLSSLKELNLVKARKQGKIVYYSLSD-"
      "-----------------";
  seqs[counter++] =
      "----------------------------------"
      "SLANLYKIFGDATRIKIVYILFEHECCVCDIAATLGTTQSNVSHQLQILKSNDLVSYRREGKQIIYSLKDSH"
      "VKDIFEKGYEHITE--";
  seqs[counter++] =
      "--------------------------------"
      "AEALAESMRAFATGSRLRLLWALLDGELTVEELAERTELSQSAVSHQLRLLRQGRLVSVRRSGRHAHYRLHD"
      "PHVVDL------------";
  seqs[counter++] =
      "-----------------------------------"
      "LADLFKVFSDATRIKILFTLLETARCVADIAEATGASQSAVSHQLRILKQAHLVTFKRCGRSIEYSLADDHV"
      "YTMLLQGMNHICE--";
  seqs[counter++] =
      "-----------------------------------"
      "LADFFTAFSDTTRLKILFELLESEKTVTEICDNTDFSQSAASHQLSKLRILKLVKVRKQGKFAYYSLDDEHI"
      "EHIIETALEHFEE--";
  seqs[counter++] =
      "-----------------------------"
      "PGVDERVAALMGALASPTRLRLFALLESGELSAGELSKAVGMRPSATSHQLRVLRDLGLVRRRREGRRCYYS"
      "LADAHLGVLLREAL-------";
  seqs[counter++] =
      "-----------------------------------"
      "LSELFKLLSDPTRMKLVLALSCGEMCVCDLGAALGMTKSAISHQLKTMKQCSVVKSRREGKNVFYSLHDQHM"
      "---------------";
  seqs[counter++] =
      "------------------------------------"
      "ADFFKVFGDPTRLKILFLLEQGEKGVNAISEELGMQQSTISQQLKLLRACRLVRFRKDGRNVLYRLNDEHIH"
      "EILALGTEHYQE--";
  seqs[counter++] =
      "-----------------------------------"
      "LAGFFSVFSDPTRLKIISALSEKELCVHELSSLLDMKQPSISQHLKMLWQARVVKKRKVGLHVFYRLDDEHI"
      "EKIYTWGYEHVKE--";
  seqs[counter++] =
      "------------------------------"
      "EAAASVATTLQALATPSRLLILTELRHGPLPVTKLAEAVGMEQSAVSHQLRLLRNLGLVTGTRSGRSIVYSL"
      "YDNHVAQLLDEAVYH-----";
  seqs[counter++] =
      "---------------------------------------"
      "FKVIADQTRMRILLALSETSLSVNEIADILTMSQSSISHQLRVLKDNRLVKGTRLGKQIHYQLTDDHIVQIF"
      "KQIIEHIEE--";
  seqs[counter++] =
      "-EHLSADPLQT-HLIT-----"
      "SGLFEPMAPQEFQTSAELFGLLSDESRLRIFWILCHYEECVINLSSLVGMSSPAVSHHLRQLKTRRLIVSHR"
      "SGKEVYYKPAD------------------";
  seqs[counter++] =
      "------------------------------------"
      "AKIFKALADDTRIKIAYILAEEELCVCDIAAIIDASTATTSHHLRLLRKMGLTKYRKEGKQVYYSLDDDHVK"
      "DLIKIAFEHQQE--";
  seqs[counter++] =
      "--------------------------------"
      "ASRLADFFALFSDCSRLRVISALAISRMCVTDLADVCRMNQTTVSHQLRSLKSMGIVESERQGKIVFYRLAD"
      "NKI---------------";
  seqs[counter++] =
      "---------------------------------------"
      "FKALGDPTRLGIVLELMETEKCVSEISSSLSISDSSTSHHLRSLRQLKLVKRRREGQKLFYSLDDHHVYLIL"
      "TIGLEH-----";
  seqs[counter++] =
      "-----------------------------------------------"
      "RMRMISLLASGELCVGNLAIALHMSDSAVSHQRKTVRALRLVGYRKQNRHVF--------------------"
      "---";
  seqs[counter++] =
      "---------------------------------------"
      "FKALGDEARLRTLEMLVNREACVSEIAEASKEQISTVSHRLKLLRAEGLVNRRREGRHIYYSLADEHVMELI"
      "HNAFEH-----";
  seqs[counter++] =
      "-----------------------------------"
      "MALLYKTFGDATRLRIMYLLLQKEMAVCDIAACLNMTHSAISHQLSVLKNLNLVKYRKFGKTVIYSLADYHV"
      "SILIATAYEHITE--";
  seqs[counter++] =
      "-------------------------------------"
      "KLFKVLSDATRLRIYALTVEEELCVCDVSASVDCSIATASHHLRTLLKQGLVKYRKEGKVVYYSLDDHHISS"
      "LVHLAMEHVNE--";
  seqs[counter++] =
      "--------------------------------"
      "ATELGEMFRLLGDPNRLRIVASCLSQPMSVGDIADTLDLSPSLTSHHLRLLRSARLLKGTRHGKQVFYDLPD"
      "CHVRQMLTNMIEHVTE--";
  seqs[counter++] =
      "------------------------------"
      "QAAAQVASTLQALATPSRLMILTQLRNGPLPVTDLAEAIGMEQSAVSHQLRVLRNLGLVVGDRAGRSIVYSL"
      "YDTHVAQLLDEAIYH-----";
  seqs[counter++] =
      "---------------------------"
      "VSEEVVIKISNFYKALSDPSRLKIMSLLNEKGLCVSCIVEKVGMTQTAVSNQLKSLRDVNLVKSERKGKNII"
      "YKLNDDHVRDILNLTMTHMEE--";
  seqs[counter++] =
      "-----------------------------------"
      "LSRMFGAFGDANRLKIMLAVADQDLCVCELGELLGMSAPAVSHHLRRLKDLSLVKTRRQGKLVYYSLDDQHI"
      "RDLLVIGQAHLQ---";
  seqs[counter++] =
      "------------------------------"
      "EETQRSVQIFKAFGDYTRYKILYLLYERELSVSEITSKIGVSQSAISHQLKLLRQTGLVSGRRDGQRILYSL"
      "ADKHIIMIFKQVKEHISE--";
  seqs[counter++] =
      "---------------------------------------"
      "FKTMSDPTRMRIILAIAQGPITVNDLAAMLDLGQSTVSHQLRLLKQARLVAGERSGKQIYYHLVDDHVLEIY"
      "ALTKAHIEE--";
  seqs[counter++] =
      "---------------------------"
      "VTPAESSRIADAFALLSDPGRVRVLALLEAGETCVCDLAEMVDMSPSALSHGLRLLRTAGVVTNRRDGRMVR"
      "YRLADSHVRLLLDVTREHL----";
  seqs[counter++] =
      "---------------------------------"
      "QRAAALFRALGDVERLKLLESLAQREVCVTELAETSRARMPTVSQRLRVLRAEGLVVQRREGKHIFYALADQ"
      "HVVELVHNALQHASE--";
  seqs[counter++] =
      "------------------------------------"
      "AKIFALLGDAGRLQLVLRCMEKPQTVGELAEASGMSQSLTSHHLRQLRDQRILASERNGRHIFYQIDDEHIS"
      "CVVRDVFAHV----";
  seqs[counter++] =
      "-----------------------------------"
      "LAETFRLLGDQSRLKILLQCMRGSVAVGDIAGSLDLSQSLVSHHLRLLRGARLVRGERQAKHIFYGIADQHV"
      "SQVLQDMAVHISE--";
  seqs[counter++] =
      "-----------------------------------"
      "LAEMFRLMGDPSRLKIIAACLGAPMCVSDIAAKYGMSQPLVSHHLRLLRAARVLRSERRGKQIFYEAADHHV"
      "KRVIGNMVEHVCE--";
  seqs[counter++] =
      "-------------------"
      "AKMNEPSILDPETLTSVSKIFKILQNEARLSIIYLLKDQELSVGEITNLINMEQSAVSHKLNALKKAHLVKT"
      "RRDGKTIFYSLDDDHVFNLLEQVITHSKE--";
  seqs[counter++] =
      "--------------------------------"
      "SEKLAEFYKTLGDKTRLRILSLLKVDERCVCELVEILGISQPTVSQHMRRLKSVHLVKERRQGRWVYYSL--"
      "------------------";
  seqs[counter++] =
      "------------------------------"
      "KTADDLAQLFSILGDGTRMRILFLLRNEETTVQSLADSLDMTHSAVSHHLRLIRPYQLVKSRKKGRNVFYSL"
      "YDDCVWHLLEEGLIHLRK--";
  seqs[counter++] =
      "--------------------------------"
      "AERLAEIMQALASPARLRILSMLSARPSTVTELSEQLQIGQTTVSNHLRLLRHLSLVTGSRAGRHIHYSLFD"
      "DHVTELLDEAIGHLE---";
  seqs[counter++] =
      "-----------------------------------"
      "LSEFFKFFGDTTRIRVIHLLLSGEMSVSAIAEKLNLEQSVVSHQLRILRTANLVKPTRDGRKIYYSLDDEHI"
      "GEIFNTGLAHI----";
  seqs[counter++] =
      "------------------------------"
      "ELLENVSDFFKALGNGTRLQIIWCLSRGELKSSELAAILQMSPSAISHQLTLLKNLKIVSVRREGKNQIYAL"
      "ADKHISQVLDSVVEHYEE--";
  seqs[counter++] =
      "---------------------"
      "TGELKHVLPEFAIETASLFGLLSDSTRLRILYLLYHREVCVRNIAEAIEMSPPAVSHHLRSLKQLGVITSRR"
      "IGKEVHYTIAD------------------";
  seqs[counter++] =
      "------------------------------"
      "EEVQGMAQMFKALSDPTRMKMAWLLDEGELCVCDMSILTKQSIATASHHLRLMKSLGIATSRKEGKNVFYSL"
      "ADHHIRTLIRMTLEHMRE--";
  seqs[counter++] =
      "---------------------------------------"
      "FSALANPIRARIVNRLTDGEASVGELSEIVGVKQPLVSHHLKVLRSAHLVSARKDGQKAMYSLIDDHVASIF"
      "LDAFNHMKE--";
  seqs[counter++] =
      "--------------------------------------FHALRSEP-"
      "RLTIMYLLLEKDMCVCELERALGMTQSAISHNLRTLRQLDLVRVRKDGRFAVYSVADEHVRALLELSRSHVM"
      "GCQ";
  seqs[counter++] =
      "-----------------------------------"
      "LTEFFKMLGNPARIRILLLLMEQDANVSDLAEQLGMTQSAVSHQLNLLKLNKLVRGCRVGKMVFYALVDEHV"
      "QMVIEKGTEHIGRC-";
  seqs[counter++] =
      "--------------------------------------"
      "YFAALSEPNRLKILYILKNGEYCPCELSEILGCTKSALSHQLRILKDKNLIKNRKDGKFIYYSIKD------"
      "------------";
  seqs[counter++] =
      "-----------------------------------"
      "LSDFFGLLSESTRLKIILELIKGEKNVSQISKNLDMSQSAVLHQLRILRQGRIVKFKKAGKNVFYSIDDEHV"
      "EGIINKAIEHL----";
  seqs[counter++] =
      "----------------------------------"
      "SVAQLFKALADENRAKIFALCQDDELCVCDVANIIGSSVATASHHLRTLHKQGIVKYRKEGKLAFYSLDDDH"
      "IKQLFTLALAH-----";
  seqs[counter++] =
      "------------------------------"
      "EVLSGLADFFSIFSDATRMKIISALTITEMCVTDISEILGINQTTVSHQLKIMRQAGVVGFRREGKILFY--"
      "--------------------";
  seqs[counter++] =
      "------------------------------"
      "EALYSEARIFKALADPNRLKIVKLLKEGELCACELTIALSSSQSTVSHHLSVLKSAGLVKERKEGKWSYFRL"
      "SEGAVIEILNQALKH-----";
  seqs[counter++] =
      "---------------------------------------"
      "FHLLRNKTRFKILLLLMEKERNVSELEEIIGGTQSAISHQLAELRNMKLVKDRQEKRMRYYSIYDSHVTSII"
      "EAAVSHINEC-";
  seqs[counter++] =
      "---------------------------"
      "ISNEVLTEIVHIHKILANPSRVKILLLLSEGQQNVSQISDAIGLEQSAVSHQLKMLKAHQLVTQARNGKAIN"
      "YQIGDSHILQLLKLSIAHAQE--";
  seqs[counter++] =
      "-------------------------------------"
      "DIFKILSEPTRLKILMALSLDSLCVCELASLLDVTQSAVSHQLRILRNAGMVDYERDGKMARYYLRDNMVV-"
      "-------------";
  seqs[counter++] =
      "-----------------------"
      "EMELINKKEITELAELFKIFSSETRLKILYLLIDTEMCVHDIAKLINMNQSAVSHQLAVLREAHLVRYERKG"
      "RVLFYSL--------------------";
  seqs[counter++] =
      "-----------------------------------------"
      "VLANSTRVQILYLLEQSELNVSELTDILKIEQSNVSHQLQRLRDYQLISQKRKGKSIYYSLDDPHIITTLNQ"
      "LMNHVQ---";
  seqs[counter++] =
      "-----------------------------------"
      "LEKIFKLLGNKQRLIILELLRERSYSVSEIINSLGMEQSAVSHQLKLLREAQLVETEKRGREVLYGLSDSHI"
      "LILLDNALKHV----";
  seqs[counter++] =
      "---------------------------------------"
      "FKLLSNPSRLQMLKVLEQRELNVGELGDLLGLEQSVVSHQLALLRKHQLVSSERVGKANYYRLDDPHILDVV"
      "NEMLEH-----";
  seqs[counter++] =
      "------------------------------"
      "EILQKMSGLLKIAGDPTRLKMLYVLVRGPKCVCDLQEEIQASQSLVSHQLKILRDNGLVKCEKIGNRALYTL"
      "SDDHVVALLSIVHEHVME--";
  seqs[counter++] =
      "------------------------------"
      "ENAETAAEALKLLAEPTRLSILALLKDNEMAVGAIAEELGRPTPAISQHLAKLRAAKLVTFRKEGTTTYYSQ"
      "KDEHVDMLVTNAL-------";
  seqs[counter++] =
      "-------------------------------"
      "VTEGMADLFKVLADDTRLKIVYALCRDELCVCDVATILGITNANASHHLRLLSHMGLAASRREGKMVFYRLQ"
      "SPHVRHLLQEVLSRGEEDR";
  seqs[counter++] =
      "-----------------------------"
      "PEYITRMSAVFQALQSDTRLKILFLLRQKEMCVCELEQALEVTQSAVSHGLRTLRQLDLVRVRREGKFTVYY"
      "IADEHV---------------";
  seqs[counter++] =
      "-----------------------------------"
      "LADCHKALGDKTRLRILALLREEDLCVGELVEILKITQPAVSQHVRKLRNARLVKERRQGQWVYYSL-----"
      "---------------";
  seqs[counter++] =
      "--------------------------------"
      "AEKIAELFKLISDGTRLRVLVLLCISEKCVSEIADAFGMSLPAVSHHLRVLKQAEIISSHRDGKEVYYSL--"
      "------------------";
  seqs[counter++] =
      "-----------------------------------"
      "LEKIFKILGNKQRLTILELLRSRSYSVSEIVDILNMEQSAVSHQLRVLREAQLVQAEKRGREVLYYLSDSHI"
      "LILLDNALKHV----";
  seqs[counter++] =
      "------------------------------"
      "ELLDDLTDLFSVFSDKTRLRIVCALAMSRTCVTDLSGVLGINQTTVSHQLRLLRNLGVVRSERDGKIIYYSI"
      "KN------------------";
  seqs[counter++] =
      "---------------------------------"
      "EQVTDIFKALSDGNRLRIMHLLIQGESSVGHIAHALDLSQSNVSHQLRILKQAHLVKGNRDGQSMIYTIDDT"
      "HVTTLLKQAIHH-----";
  seqs[counter++] =
      "-----------------------------"
      "PKIFEMSTDLFSILSDVSRVRILWLLCHTEDCVANIADAVDMSSPAVSHHLKLLKSANILKYTKKGKEVYYT"
      "LAD------------------";
  seqs[counter++] =
      "------------------------------------"
      "AAMFKLLGDPTRARLLALLEAGELCVCDLAAATGTQEATVSQSLRMLRASGVVTGRRQGRLVFYRLADAHV-"
      "--------------";
  seqs[counter++] =
      "---------------------------------"
      "ERLAGRFKGLADANRIKIAYLLTREELCVHDIARLVGISIANASHHLRLMRSLGVTKTRKKGTTVFYSLADR"
      "HVHTIVLLGMEHMEE--";
  seqs[counter++] =
      "------------------------------"
      "QVIVNLSSLYKVFADKTRLEILYALHENEMCVCDLAVLLNMTKSAISHQLKTLRLANLVKNRRVGKVVYYSL"
      "ADEKVYEIFNQSFKQLTE--";
  seqs[counter++] =
      "---------------------------------"
      "QKVSQLYKVLSDPTRLKILLYLKQGELNVTALSEKLNMEQSAVSHQLKLLRENHVVKTDRVGKTIFYILDDH"
      "HVLDILNQTIQHI----";
  seqs[counter++] =
      "----------------------------------"
      "ALSDFFRIFGDQTRLRILYALAKTELCVCDLAKLLGASQSAVSHQLQVLRSHRLV-----------------"
      "----------------";
  seqs[counter++] =
      "------------------------------"
      "DVAEKLAGAFKLLSVEARIRIVQVLKRRALCVTELTSQLGISQSATSQHLRVLKDARIVKFQKRGFHVYYHL"
      "--------------------";
  seqs[counter++] =
      "--------------------------------------------"
      "DPTRLKILLSLKEGELCSCDISEISKISISATSHQLRLLRDRKLVKYRKEGKFVYYELYDEKI---------"
      "------";
  seqs[counter++] =
      "---------------------------"
      "ISDHTATHLADTFSLLGDPSRIRVLGTLLDGPKRVLDIAQACGHTQSATSHSLRLLKAHHVVAGERHGREIH"
      "YALADDHVRALLTLALAHI----";
  seqs[counter++] =
      "---------------------------------"
      "QKILEILKILSDETRLKIVSLLAENELCVCELMEALRMSQSRISNHLRILRNTRIIEAKREGKWIFYSL---"
      "-----------------";
  seqs[counter++] =
      "-----------------------------------"
      "LADFFDIFGDTSRIKILLALHDKSLPVSSIAELTGLSASAVSHHLSLLRGRRVVKVERKGKYRVYELDDDHV"
      "SSVLKMAISHIQEVK";
  seqs[counter++] =
      "---------------------------------"
      "EGIAVIFKALADDTRLKIVYALSQAELCVCDVAALINSTKSTASYHLRLLNHMGLAKFRKDGKLVYYRLADQ"
      "HI---------------";
  seqs[counter++] =
      "--------------------------"
      "AIDAEAVQGVSALFKALGDETRLKVLALYKGEELCVCDVANIVGSTVATASHHLRLLRNIGIANYRKEGKLA"
      "FYSLRDAHI---------------";
  seqs[counter++] =
      "------------------------------------------"
      "LGDPTRLKIIYTLSETSMCVSDIAKTLDLSQSLVSHQLALLREAELVKVKRVSRNAIYSLDDAHVLTIFKQA"
      "HEH-----";
  seqs[counter++] =
      "-----------------------------PQLAAA-"
      "ANTFAMLASPARLHLVTLMSGGRFDVGTLAEKVGLSLPTTSQHLSKLRLAGIVSARRAGRHSYYTVEDPHVL"
      "SLVEQIFEHI----";
  seqs[counter++] =
      "---------------------------------------"
      "FKVLGDPTRTKIVLALDNREVCVCTLADTLGMTKSAVSHQLAILKANNIVKSRRDGKQVYYSFDDEHITDII"
      "EIAQAHIKD--";
  seqs[counter++] =
      "------------------------------------"
      "ARLFKGLGDETRFRIYALYLEAELCVCDVAGILGTSVATASHHLRLMKNLGLTRSRKEGKMVYYSLDDDHVR"
      "LLVKLAIDHAAE--";
  seqs[counter++] =
      "---------------------------------------"
      "FQALGDPERLRLMIRLSEQEICVSELAELAQEQLTTVSARLKSLYAARLVKRRRQAKHVFYSIADDHVLQMI"
      "RGAVAHAAE--";
  seqs[counter++] =
      "-------------------------------------"
      "DIFKIIGDPTRLMILHAIEFHELCVCDLGHLVGVTKSAISHQMKLLKKYGLVKGRKVGKMVYYSVIDDNVKN"
      "LIHAGYNHV----";
  seqs[counter++] =
      "-------------------------------------"
      "ETFRLLGDPTRLKILLACLSEPKCVNDVASEVGITGSLTSHHLRLLRGARLVRAERQGRQIFYVAADSHVNA"
      "MLAEMVAHIRQ--";
  seqs[counter++] =
      "-----------------------------------"
      "LAEVFKVFGDSTRIKILYDLFEGEKNVTEICQDLEMNQSAISHQLKILRTARLVSGKRMGKSILYSLADEHV"
      "KTIIAMGIEHIEE--";
  seqs[counter++] =
      "---------------------------------------"
      "FSALADRSRLKILYALSETELCVCDVASLLGMKIATASHHLRKLRDLQILKYRNDGKLAYYSLKDQRVAEIL"
      "HHTLNQLVE--";
  seqs[counter++] =
      "-----------------------------------"
      "LADLFKVFSDSTRMKIMYKLFDGEVSVGQIATSLDMSQSAISHQLKYLKESNLVKSKRNGKSMLYYLADDHV"
      "KIIIKTGLEHIEE--";
  seqs[counter++] =
      "-----------------------------------LVEFFKTLGDFTRLRIV-"
      "LELKTKRCVGELAEELEMSHSAVSHQLNILKANGIVKSQRQGKYIYYIVQDEYV----QNAIE------";
  seqs[counter++] =
      "----------------------------------"
      "AVAEVFKLISDGSRLRILWLLCHREVSVGDIAEMMDMSNPAVSHHLKLLKQSGLVDSRREGKEVFYRLAD--"
      "----------------";
  seqs[counter++] =
      "--------------------------------------"
      "FFKQLGYSTRVRILCYLIQEPRKVSDIAVHLNMTLSSVSHQLRVLREAGLVSGQRQGKTITYQKKDDHVATI"
      "IQNTLDHL----";
  seqs[counter++] =
      "---------------------------"
      "VSPDVLELIAERFRVLADPARLQILNVLQGGEQTVTELMRTTGFRQAKVSKHLQLLYNLGFVDRRKEGLHVY"
      "YRL--------------------";
  seqs[counter++] =
      "---------------------------------------"
      "FKMLSDKTRLSIMLLIKEQEMNVSEISRALNMEQSAISHQLSALRSERLVKSRREKRSVFYSPNDQHVYDIL"
      "TQVIDHLETC-";
  seqs[counter++] =
      "------------------MARAGDV----"
      "PSDYEPVSALFKALANPVRAAIVHLLSDRERTVGQLVEALGLPQPLVSQHLRVLRGALMVATRRQGQEIWYS"
      "VCDQHVAHILGDAMKHTQE--";
  seqs[counter++] =
      "---------------------------------"
      "KDIAEFFKLFSNDGRLKIISSLATDNLTVNEIVERTKLSQSLVSQQLKLLKNARILTNEKIGKTVTYSIYDR"
      "HILHLLKDVAEHLDE--";
  seqs[counter++] =
      "----------------------------------"
      "SVAQMLKAIADENRAKIYALCQDEELCVCDIANIIGITVANASHHLRTLHKQGIVRYRKEGKLAFYSLDDEH"
      "IRQIMMIVLEHKKE--";
  seqs[counter++] =
      "-----------------------------PDVGE-"
      "MVQIFKALADETRLRIYSLTLESEMCVCDVAAVIQSSSATASHHLRYLREHSLAKSERRGKMVYYALADKHV"
      "ADLYEHAIEHTME--";
  seqs[counter++] =
      "---------------------------"
      "VSQEVVQQVAEYFSLLSEPMRLRLLHLLRDEEKCVQELVDATQTSQANVSKHLKVMWQAGILSRRSEGTSAY"
      "YRVEDEMIFELCNRVCDRL----";
  seqs[counter++] =
      "---------------------------------"
      "EKVSQLFKMLSDPTRLKILLYLKDGEQNVTAITQAVEMEQSAVSHQLRLLRENHVVKSHREGKAILYSLDDS"
      "HVLDILNQTLKHVEQ--";
  seqs[counter++] =
      "-------------------------------"
      "IAEHAAEVLKAIAHPVRLQIVELLQAEEMCVGDIVNALGAKQAITSQQLNMMKAKGVLSCRRDGARVYYRIE"
      "NRNVIKLLDCIYDHCEK--";
  seqs[counter++] =
      "---------------------------------"
      "ERIAETLKALSDPTRLRIVSLLRHGELCVCDLTEALQTPQSKVSRHLAFLKNAGWVRARRSGKWVYYQILD-"
      "-----------------";
  seqs[counter++] =
      "-----------------------------------"
      "LSELFWVLSDATRIRILYALSEKEMCVCELARLLNNKQSSISHKLRILRNSKLVGFKRNGK-----------"
      "---------------";
  seqs[counter++] =
      "--------------------------------------------"
      "NPTRLRILLLLSKKDMCVGKIAEILRMDQSAVSAQLKVLRHLNLVKAKRHGRYMRYKLNNKHV---------"
      "------";
  seqs[counter++] =
      "------------------------------"
      "EHSQIAAETFRMLADATRVRILWALFHDELSVNALAEHVGAVPTAVSQHLAKLRLAGLVSSRREGTFVYYSA"
      "SDAHVKALVAQAL-------";
  seqs[counter++] =
      "----------------------"
      "GEIMVILPDRLQAQSELFKPLSDPTRLKILYLLRNGELCVCEIIFALKKPQSTISHHLNILKKAGFIKGRKE"
      "GVWIHYRLADAEIVGVIDNLTSILNE--";
  seqs[counter++] =
      "---------------------------------------"
      "FKALADETRVKIIALLQEPNLCVCDLAQITELTISGASHHLRLLKNMGLARSHKEGKHVRYCIHDEHVKIIL"
      "EEALNH-----";
  seqs[counter++] =
      "------------------------------"
      "EQLQSLADTFKLMGDKTRLTILALLRERELCVCDLVDVTGMSQPSISQHLRKMKDAGLVSETRKGQWIYYSL"
      "--------------------";
  seqs[counter++] =
      "-------------------------------------------------------------------"
      "ALGTSESAVSHQLRRLRDQNLVLPRKEGRVVYYRLADAHVTDLLRNVLEHVGE--";
  seqs[counter++] =
      "-------------------------------------"
      "DLLGALANANRLKILSLIIDGELCVSAINAHVDLSEAALSHRLAKLRKLRLIESRRQGTTIYY---------"
      "-------------";
  seqs[counter++] =
      "------------------------------"
      "EQAEELAGMFHLLGDVNRLRLICACLEEAVCVQDLADRFSLTPSLVSHHLRLLKAARLMRAERRGKQVFYTV"
      "NDEHV---------------";
  seqs[counter++] =
      "------------------------------------"
      "AKIFALLGDAGRLQLVLRCMEKPQTVSELAEATNMSQSLCSHHLRHLRDQRILASERRGRHILYRIDDDHIS"
      "RVVRDTFAHVHE--";
  seqs[counter++] =
      "-------------------------------------"
      "EVFRMLADSTRIQLLWALIDRELSVNELAGEVGKLPASVSQHLAKLRMSRLVHTRREGTQIYYRLENEHIAR"
      "LVTDALD------";
  seqs[counter++] =
      "---------------------------------"
      "KKIAEFYKALGDEVRLKILQMLSEQEMCVCEIIERLDMSQPAVSHHLKILRQVGLVKDSREGKWIYYSLHD-"
      "-----------------";
  seqs[counter++] =
      "-----------------------"
      "ELPALDDQHIASLAHLFHLLGDEGRIRLVLACMAGPVPVSELSAVTGMSQSLTSHHLRHLREARILRSERQG"
      "KQILYRLDDHHI---------------";
  seqs[counter++] =
      "-------------------------------------"
      "DIFKLLSHPMRLQIIYMLEQQTMNVGEIVERLGLEQSAVSHQLTLLRKGHLISTCQIGKIVCYSLNDKHILD"
      "IVNEALEHTQ---";
  seqs[counter++] =
      "------------------------------------------"
      "LSDPTRLRLLALLVREELSVAELQEILGMGQSRVSSQLALLRQVDLVTDRRDGKKAFYSIRSNRTLALLKSA"
      "IDSVSE--";
  seqs[counter++] =
      "------------------------"
      "IQPLSRKNAQSAEKLFKCLASASRLKILFVLLESEKSVGDITVDCDMSQPLVSQHLRHLRDNNLVYTKRHGK"
      "QVYYSIADEHIKHVVADCIQHVQ---";
  seqs[counter++] =
      "---------------------------IASELAPAVA-"
      "LFRSLGDPARLAILDRLARGEARVVDLTDELGLAQSTVSKHLACLRDCRLVDFRVEGRQSFYALARPELPAL"
      "FRSA--------";
  seqs[counter++] =
      "------------------------------"
      "EEAAQLTSVLSLMADPTRARVLALDMVKELCVGDLALALESNEDAVGYALRLLRTAGLVTNRKQGRVVFYRL"
      "--------------------";
  seqs[counter++] =
      "------------------------------------"
      "ADLLKALSNPGRLRILCALVPGEMSVGDLETALGASQSYVSGQLLRLRNEGLVSCTRDGRSIRYQLAD----"
      "--------------";
  seqs[counter++] =
      "------------------------------------------"
      "LSTPSRLLILARLREGPLPATELAAEVGMEQSACSHQLRLLRNLGLVVGERRGRSVVYALHDDHVAGLLDQA"
      "VYHVE---";
  seqs[counter++] =
      "------------------------------------"
      "AEYFKAISDPARVRIIYALANGELCVCELMLIMGMQQTVVSHHLKILKYANIVSDRKSGKWVNYSLADRRVL"
      "--------------";
  seqs[counter++] =
      "------------------------------------"
      "ANIFKLLGDYNRIRIINALKIKELCVCELSILLDMSQSSISHQLRILRHHNIVKNRKENKRVFYSLNNDKIF"
      "KLIEESI-------";
  seqs[counter++] =
      "-------------------------------------"
      "ELFTQLSSSTRIRLLCILAIEEMCVCELADMLKMSQPSISHHLRLLRQSGVVKYKKSGKRVIYYITD-----"
      "-------------";
  seqs[counter++] =
      "----------------------------------"
      "TLAELFKALGDPTRLNVLQLLTERQLCVGAIARRLGVTQPAVSQHLKVLKHLGLVKASRDGYHIHYSI-"
      "NQDMLASYKTHIDEWQ---";
  seqs[counter++] =
      "-----------------------------------"
      "MADILKVLGDPNRLHILSLISRQELCVCEITSILNISQSNASQHLARLRSVDLVKERRNAQWIYYSL-----"
      "---------------";
  seqs[counter++] =
      "---------------------------------------"
      "FQALSEPIRLQILDLLQEQELCVCEIREKIKISQSKLSFHLRILREAKLARSRQQGRWVYYSLNPEQLL-"
      "LLEQYLNQLRE--";
  seqs[counter++] =
      "-------------------------"
      "KALQREHLQSAADLAKSLSDENRLRILDCISRGKQSVGGIAKELSLSQPLVSHHLRELRRTLLVKIERNGAF"
      "VYYELSD------------------";
  seqs[counter++] =
      "-----------------------------------------"
      "ILGEPSRLKIVLALSEGDMCVYHIVKAVNSNQSAVSHQLRILRENKVIKSFRKGQNIVYSLDDEHIMQIINI"
      "VKTHVEE--";
  seqs[counter++] =
      "----------------------------------"
      "SLAKYGKAISDPKRIELMDLLVQAEKNVDVLSKETGMSIASTSHHLQILKEARLVSDRRKGRNIFYQIED--"
      "----------------";
  seqs[counter++] =
      "------------------------------"
      "EHAEDLAQVMQALSSPGRLLILARLDDSPCSVSTLVEDCGMAQATISNHLRILRHLDLVTGQREGRQVIYSL"
      "YDAHVQEFFRQALGHI----";
  seqs[counter++] =
      "-----------------------"
      "DLSSLSRSQAEVASELFKSLSNPNRLQIVAALALGEHPVGDLETMLGIKQPTLSQQLAELRDAGFVESRREA"
      "KQVFYRLGDKRLLAL------------";
  seqs[counter++] =
      "-----------------------------------"
      "LADFYKLFSDSSRIKILFVLLSGAHCVKHIAEKAEMSQSAVSHQLAVLRRSNIIRQTRSGQNITYSLADDHV"
      "KLLLELAIAHIRE--";
  seqs[counter++] =
      "------------------------------------"
      "ANLFKALADETRLSIYALTIEEEMCVCDIAAVIGSSMATASHHLRYLRERSLAKSERKGKQIYYSLSDNHVR"
      "QLVKIAHEHTKE--";
  seqs[counter++] =
      "------------------------------"
      "QAATQTASLLKTLGNPDRLLLLCQLTQGEACVSDLEASLGIQQPTLSQQLTVLRNEELVATRREGKRIYYSI"
      "AD------------------";
  seqs[counter++] =
      "------------------------"
      "LQMNMTQAATQTASLLKTLGNPDRLLLLCQLTQGEACVSDLEASLGIVQPTLSQQLTVLRNEGLVATRREGK"
      "RIYYSIADEKLFTL------------";
  seqs[counter++] =
      "-----------------------"
      "EMQDIAAQLQELHARVCKAIADPKRLLIINELRDGELSVGDLCEALGFSQSNASQHLSVLRERGIVNARRSG"
      "NNVFYSLRSRKIV----QAVDLLRE--";
  seqs[counter++] =
      "------------------------------------"
      "AKFFHGLANPTRLKIVETLLAGEMSVSQIVDAVGVSQSQVSNQLACLKWCGYVTSRKEGKYILYRISDERVR"
      "AILQLA--------";
  seqs[counter++] =
      "------------------------------------"
      "AEGFRLLADPTRIKILWALLQGESSVACLAEMVGAAPTAVSQHLAKLRLAGLVKGRREGTYVHYSAADGHVR"
      "ALLAEALFH-----";
  seqs[counter++] =
      "---------------------------------------"
      "FKALSDFNRVRIMEFLENGEASVGHISHSLNMTQSNVSHQLKLLKSTHLVKSKRQGQSMIYSIDDIHVSTLL"
      "KQAIHH-----";
  seqs[counter++] =
      "----------------------------"
      "APEVGVTLARALLALTDPSRVRLCRLIARQAMTTADLADRLTMTRPQVSRHLRALRELGLVRMERHGRHVLY"
      "EL--------------------";
  seqs[counter++] =
      "---------------------------------"
      "QMLDEFFKSLSEPVRLRVMYLLERGELCVCDIVSSLEVSQSVVSRHLAYLRNAGLVSSRRQGVWIYYQL---"
      "-----------------";
  seqs[counter++] =
      "------------------------------"
      "DLAASVSEKLKVYAQPQRLMILSCLWRGERTVADIGQATGIVQPALSQQLAELRRADLVQTRKEAKQVWYRL"
      "AD------------------";
  seqs[counter++] =
      "---------------------------------"
      "KKLASFFDVLSDGTRLKILSALAITPMCVSDLSAVLEINQTTVSHQLARMRLAAMVDFRREGK---------"
      "-----------------";
  seqs[counter++] =
      "-----------------------------------"
      "LAELYKLLGNVTRLKILLALAQGELCVCDVAHVLGLTVAATSHQLKLLRDQGWLAMRNDGKMVYYRL-----"
      "---------------";
  seqs[counter++] =
      "-------------------------------------"
      "KILSLLKNPVRLQILYILSQQSLSVSEIVELLHLDQSLVSHHLSDLRKYQLVSTKRDGKSIFYELDDPHILD"
      "IVNETLEH-----";
  seqs[counter++] =
      "-----------------------------------"
      "LSETFRLLGDPSRLRILLHCEEGPKSVTDISETLELSQSLVSHHLRLLRGARLVTRVRHSKQMFYEISDQHV"
      "GDVLLDMLSHVRE--";
  seqs[counter++] =
      "-----------------------------------"
      "LAELFAALSDPTRLRLLNLMRDREVCVCDFVEILGQSQPKISRHLAYLRRAGIVCARREGKWMHYRIE----"
      "---------------";
  seqs[counter++] =
      "---------------------------------------"
      "FQQLGDPTRLKILWILCHCRECVSDIAAAVGMSDAAVSHHLQLLKRSGLIVGSRVGKEIHYTLSDERRAGLL"
      "HRMMDALFE--";
  seqs[counter++] =
      "-----------------------------------"
      "MADIFSVVADPTRRELLGTLLSAELSVGQLVERLGVSQPTVSKHLRVLRDIGLVTSREEGQHRYYRL-----"
      "---------------";
  seqs[counter++] =
      "------------------------------"
      "EYFQTVALVFRQLSDANRVRLFWLLCHCEECVVNLAAMMGMSPPALSHHLRQLRESGLIVSRRDGKEVYY--"
      "--------------------";
  seqs[counter++] =
      "------------------------------------"
      "AEVFNQLSDGTRLRILWLLCHSEECVNDIAAAVRMTAPAVSHHLKTLKQNGIIKSRRLGKEVLYTLED----"
      "--------------";
  seqs[counter++] =
      "-------------------------------------"
      "QIYDALSDFTRFQILGALLLGEKSVTELQELLSVSQSATSHQLRLLRDRGLVTAKRDGRRVIYSLADDHVIT"
      "LISVGLAH-----";
  seqs[counter++] =
      "--------------------------------------"
      "FLNLIADKRRIDIIYLIMKKRLCVQDIAEIINETVANTSYHLQQLKKGNIVKVEKEGKEVFYSLSDKHVYEI"
      "LENVLEHI----";
  seqs[counter++] =
      "-----------------------------"
      "PQLLETAAGTLRMLAEPTRLNLLFQLTDGPKTVTELTAAVDVPRTVVSQHLAKLRLSGLVDTRKDGRHVIYS"
      "LHDGHLIRLIRETINH-----";
  seqs[counter++] =
      "--------------------------------"
      "AQAATALLKVLANENRLMILCTLMGGEMSVGELNTAVPLSQSALSQHLASLREAGLVSTRKEAQTVYYRLQ-"
      "------------------";
  seqs[counter++] =
      "------------------------------"
      "QAADAAVELLKALANPVRLKLLCFLVEQERSVGEIASRLGVRETLVSQHLSLLRRDKLVAYRRDGQTLWYRL"
      "AD------------------";
  seqs[counter++] =
      "-----------------------------"
      "PKTSELQAKLFRGFADPSRLAILETLRDGPLTVGEIVQATGLSQSNVSNHLGCLRDCGLVTATQQGRFVSYA"
      "LSD------------------";
  seqs[counter++] =
      "------------------MAESQEMH---"
      "QEEAQNAASFLRSLGNPHRLQILCRLALGEQSVGQLHQFFDLSPSAFSQHLAVLRQQQLVSIRKESQTVYYS"
      "IKD------------------";
  seqs[counter++] =
      "------------------------------------------"
      "LSDAGRLRLLLWLAQREMCVSELVALEQDKVSSVSARLQMLHAVNLVTRRREAKHMFYALADVHVHRLLRNI"
      "LDHAAE--";
  seqs[counter++] =
      "---------------------------------------"
      "FAVIAEPSRRRILDLLLQSESNVTDLAQALGLSQPLVSKHLRTLRQSGLVKVRKSQQHIY------------"
      "-----------";
  seqs[counter++] =
      "------------------------------------"
      "ARFFRVLGDPVRMGILELLLEGEKNVSEIVSRLGMSQSRVSNHLACLRWCGLVSVRRKGSFIYYSLADEQLR"
      "ELLEIANDRVEK--";
  seqs[counter++] =
      "----------------------------------"
      "TLSDILHLMGEPNRLRLLVTCLEGAKSVSELAQQLQLSVPLTSHHLSLLRSARLLVANREGKHIYYSIYDAH"
      "VRCILDDMLKHFTE--";
  seqs[counter++] =
      "---------------------------------------"
      "FELLSDGNRLRLLLCLHHADICVGDLAAALDMTGTAVSHALRLLRNQGWVSATRDGRSMRYRLTD-------"
      "-----------";
  seqs[counter++] =
      "-------------------------------------"
      "DFFKALAHPTRIVLVEDLAEGEKCVCDLAQKIDADISTVSRHLRELRNAGIVANQKRGNQVFYSL-------"
      "-------------";
  seqs[counter++] =
      "-----------------------------"
      "PEVDREIIEFLKALSNPIRLKILKLTRDNWLCVCLLSEVLGEDQTLISHHLRTLKTLDLVKERREGRMRFY-"
      "---------------------";
  seqs[counter++] =
      "------------------------------"
      "ENADQAADFLSALANNKRLLILCKLLHNEMSVGALAKAIDLSQSALSQHLAKLRALDLVSTRRDAQTIYYMV"
      "SSPHI---------------";
  seqs[counter++] =
      "---------------------------------"
      "KNLSSYFKGLADENRLRILNLLFHGELCGCDIQYVLGASQSNVSRHLSYLKNAGLVNDRRKANRVYFSL---"
      "-----------------";
  seqs[counter++] =
      "------------------------------------"
      "AELFKVLSSATRLRLLRTLAEEVSTVSRLAERSGLAQPLVSQHLRTLRSAGLVSVERVGREAHYSVADTHVT"
      "HIVEDAVHH-----";
  seqs[counter++] =
      "-----------------------------------"
      "LAETFGILSDSTRLSIVLACMETEVSAGDIATKLKVSPSLVSHHLRLLRAVRIVRSERRGKQVFYTMTDACV"
      "RDILTTMINHLPE--";
  seqs[counter++] =
      "------------------------------------"
      "AKFFRGLADPSRLALLLALRPGEKTVSTLSEETGLSQSNVSNHLACLKDCGLVVNRQEWRHVYYRIADSKVL"
      "TL------------";
  seqs[counter++] =
      "---------------------------------------"
      "FKCLSDPTRLQILHLLMEGEHCVGDIALKIGTTQANISKHLSLLKNAGLVVSNKQGMKVIYSLQ--------"
      "-----------";
  seqs[counter++] =
      "----------------"
      "VPLPTSTLMQDIAAQLQELHARVCKAIADPKRLLIINELRDGELSVGDLCEALGFSQSNASQHLSVLRERGI"
      "VNARRSGNNVFYSLRSRKIV----QAVDLLRE--";
  seqs[counter++] =
      "---------------------------------"
      "KSLNEQFSVIANPNRLEMLEFLAQCEYSVDDLAKVMGLSVANTSHHLQQLRLAGLVASRKEAQRVFYRLKGD"
      "GVVEL------------";
  seqs[counter++] =
      "------------------------------------"
      "AEMLKALADPHRLGILLRLSKRELSVGELAEIEQEKVTTMSARLKVLLTAHLVKRRRKGQSVLYSLADTHVL"
      "NLVDNAIEH-----";
  seqs[counter++] =
      "-----------------------------PDVVQ----"
      "LFKVLADETRLEILRLLALTDLRVGEIVAHLGLPQNAVSYHLKQLRRLRLLRDHRSARDIYYSVDLDHLQAL"
      "YAAA--------";
  seqs[counter++] =
      "-------------------------------------"
      "ETLKTLANQKRLEIVQLLGQGELTVSEMIEMLGISQSNLSQHLAVLRRYQIVATRKEGLYVYYRLTDGHI--"
      "-------------";
  seqs[counter++] =
      "---------------------------------------"
      "FKILGDENRLRILNLLRKGELCVCEIELVLETTQSNVSRHLGKLRNEKIISFEKKAQWIYYRI---------"
      "-----------";
  seqs[counter++] =
      "------------------------------------"
      "ADLFKALSSPARLRILSALIAGPSDVGSLADATELSQPLVSQHLRTLRLAGIVQVERIGRNAVYSLHDEHIA"
      "HIVGDAVSHVSE--";
  seqs[counter++] =
      "------------------------------------"
      "ARLFKVLGSESRLALLRILQAKPATVGVLVEKSGLTQPLVSQHLRVLRQTGLVTRDRQGKEVTYQIADHHVA"
      "HLIDDAIIH-----";
  seqs[counter++] =
      "-----------------------------------"
      "LANLFKVFSDSTRIRILFSLFDYEKNVNTISKELNLSQSAISHQLRYLKDSNLVKSQRDGQAMIYALSDRHV"
      "KFIIKLGLEHLYE--";
  seqs[counter++] =
      "-----------------------"
      "KVEAVKEELAASISQLFKVLADERRFKILYALTKQELCVCDVALIIGATVATTSHHLRTLSKQSILTHEKIG"
      "KMVYYQLSNPMIQQLVLDAMNQEKE--";
  seqs[counter++] =
      "---------------------------------"
      "KAVVRIFDVLGNRTRLRILLALASEELCVCDIAHALNLSISAASHQLRALHDRDWLRMRNDGKMVYY-----"
      "-----------------";
  seqs[counter++] =
      "---------------------------------"
      "ENVGEFAKALGHEKRLLIIELLSSHERCVEDLATAMGIGVKSVSAHLKVMRTQGILTTRKEGLRVYYRLRND"
      "NILKLFQ----------";
  seqs[counter++] =
      "---------------------------------"
      "KELSNFFNAFGNPTRLKILLALKEEELCTCDLSNITGLSVSAISHQLRVLKGRKNVNYRGDGK---------"
      "-----------------";
  seqs[counter++] =
      "---------------------------------------"
      "FHALGEPLRLKVIEILHREELCVCDLCERLHLRPSKLSFHLRALRQANLVLSRQQGRWVYYRL---------"
      "-----------";
  seqs[counter++] =
      "-------------------------------------"
      "EVLRLLADRTRLAILAMLDGTEMPVNAIAEALGRPAPAVSQHLARLRAGRLVTSRRDGTTVFYGQPDEHVAA"
      "LVANVLQHTEPHR";
  seqs[counter++] =
      "-------------------------------------"
      "EVFKAVADPCRLRIVKLLKEGELCVCEIMTALDKPQSTTSHHLSILREAGLVRERKDGKWSYYRLAD-----"
      "-------------";
  seqs[counter++] =
      "---------------------------------------"
      "FRALSDPIRLNVINLLQEKEMCVGDICLALKIAQPKLSFHLRVLRESGLLQTRQEGRWIYYRL---------"
      "-----------";
  seqs[counter++] =
      "---------------------------------------"
      "FKALADANRRKILFLLKESDLTAGEIASEFDISKPSISHHLNILKNAGLVEARREGQQIYYSL---------"
      "-----------";
  seqs[counter++] =
      "------------------------------------------"
      "LGEVNRLSLLALLHAGDLCVSDLAVAVGMSDSAVSHALRLLRAHGMVTAHREGRLVRYRL------------"
      "--------";
  seqs[counter++] =
      "---------------------------------------"
      "FKILSDETRLRIIILLAQEELCVCQISGVLNVSQPKVSKSLSRLRDLNLVIDERKEKFVFYKLKTENFVSTI"
      "RNIMDNLNESR";
  seqs[counter++] =
      "----------------------------------"
      "TLESLFSALADRTRLEIVLFIMRGKASVQEIARGINKSQSLVSHHLACLRNCGIVKTERKGKYVYYSLLDNE"
      "VVSIIKLAVEH-----";
  seqs[counter++] =
      "---------------------------------------"
      "FQALADPSRRAIFESLTRGEAAVKDLTTRFDISQPAVSQHLAALKDAGLVSGRREGRHVYYRVE--------"
      "-----------";
  seqs[counter++] =
      "------------------------------------"
      "AALFHALSDAGRLRTLAILAEQSSSVSHLAEVTGERIGTVSARLKVLLQANLVTRRREGQSAIYSIADQHVL"
      "ELIHNALEHVNE--";
  seqs[counter++] =
      "-------------------------------------"
      "DFFGIVADETRLRIIGLLNQKELCVCEMCEILGLSQPKVSRHLSKLRDAGIVIDSRQGQWVFYYL-------"
      "-------------";
  seqs[counter++] =
      "---------------------------------------"
      "FKALGDNNRLRILSMLNVRELCVCEINAVLKVSMSTISSHLKILRNAGLVTSRKDGRWIIYRLE--------"
      "-----------";
  seqs[counter++] =
      "------------------------------------"
      "AELFKALATPSRLKILLTLSHGPASVSNIVIATELSQPLVSQHLKVLRGIHLVSVQRDGREAIYSLMDDHVA"
      "HIILDAMAHVNE--";
  seqs[counter++] =
      "------------------------"
      "MMSVSEETSDEAARQLKAVADPVRLRILYALSKEPLCVCELSVLLNMSMPAVSHHLRILLSAGLLKVRKEGK"
      "FACYHLRDSH----------------";
  seqs[counter++] =
      "-----------------------------"
      "PQAAHEASDLLKALAHHTRLLILCILAKQERTVGEIENILGIQQAMVSQQLARLRLEGLVNTRRQGRLVYYS"
      "IGNVSVLAFLESLFD------";
  seqs[counter++] =
      "------------------------------------------"
      "LAEPTRLHLLWQLSNGPKTVTELTDASGAARTVVSQHLAKLRLSGLVDTRKDGRHVIYSLHDGHLVRLIRET"
      "INH-----";
  seqs[counter++] =
      "----------------------------------"
      "TISQIFKILSDETRVKIVALLTENELCVCDLANIVEATVAATSHHLRFLKKQGIANYRKDGKLVYYSL----"
      "----------------";
  seqs[counter++] =
      "------------------------------------------"
      "LGDPTRLRVTALLSGEELCVCDLAWVVGLAQNLVSHHLRLLKGAGLVTGRRHGRLVMYAL------------"
      "--------";
  seqs[counter++] =
      "------------------------------------"
      "AAIFKALGEVNRTRIVKALSLEELCVCDIATIIDATIATTSHHLRSLHGQGIIKSRKEGKMVYYSLDDDHIR"
      "QIVSMAFLHQEE--";
  seqs[counter++] =
      "------------------------------"
      "EDAEIAAGFLSAMANPKRLLILDSLVKEEMAVGALAHKVGLSQSALSQHLSKLRAQNLVSTRRDAQTIYY--"
      "--------------------";
  seqs[counter++] =
      "------------------------------------------"
      "LSDPHRQSILKMLAHQEMGACEIIHSIGLSQPAVSHHLKILRQARLITSQKQGKMVFYSL------------"
      "--------";
  seqs[counter++] =
      "------------------------------------"
      "AEYFKALSHPTRIKIIELLSKKEMCVCQMMAALNLDQSHVSRHLMVLRANEMVKTRREGTIIFYSLTDENII"
      "--------------";
  seqs[counter++] =
      "-------------------------------------"
      "EMFRAFSDRTRLRILNLLLRGEMCVGDLVSILEMSQPRVSQHLSCLRNSGLVVGRREGQWNHYSL-------"
      "-------------";
  seqs[counter++] =
      "---------------------------------------"
      "FKLLSNPTRLNILMLLEHEQLSVNEIVTQLEITQPQVSHQLAILKEQQLVSAKKIGKKSLYQLSDPHILSV-"
      "-----------";
  seqs[counter++] =
      "-------------------------------------"
      "DIFEALSDPHRRKILDMLKHGELCSSDIASQLDITPASVTHHLNKLRSANLIIKTRKGRNIYY---------"
      "-------------";
  seqs[counter++] =
      "------------------------LERLEPQISEA-"
      "ARLMEMLSHPARLRILCTMLGGEKSVQELAINASLSQPAMSHHLRKLRDSELVNTRRDKQTIYYSLKGEHVA"
      "AVLE-VLEHL----";
  seqs[counter++] =
      "---------------------------------"
      "QARAERLRALGEPTRLRIYALHAGAELCVCDLAWIIGSSQGLVSHHLRQLRAAGLVTSRRDGKLVMYRL---"
      "-----------------";
  seqs[counter++] =
      "------------------------------------"
      "ARIFKVLGDRNRTAIYALCENDTLCVCDIATIIDASVATTSHHLRTLYKEGVVTYEKKGKLAMYALDDNHIR"
      "QLMMTTLEHAEE--";
  seqs[counter++] =
      "---------------------------------------"
      "FELLSDANRLRLLLCLHHAPICVTDLSVALGMSGTAVSHALRLLRSQGWVSATRDGRSMRYQLAD-------"
      "-----------";
  seqs[counter++] =
      "------------------------------------------"
      "MADPLRLQVLNLLSKQELCVCDLCDRLQVKQPKLSFHLRQLREAGLIQARPQGRWTYYSL------------"
      "--------";
  seqs[counter++] =
      "-----------------------------"
      "PDFVETSAALLQAMANPARINILIILAEREVSVGPLSELVGLSQSALSQHLAKLRQAGLVSTRREAHTVYY-"
      "---------------------";
  seqs[counter++] =
      "------------------------------------"
      "ADLFKVLSNPVRIQILDALRLGEQSVGYIAEWLEIEASAVSQQLAVLRSRNLVTSRKQGNYVFYSVRD----"
      "--------------";
  seqs[counter++] =
      "---------------------------------------"
      "FQALSDPLRLQILQLLRHQELCVCELRDHLDIAQSKLSFHLKTLKEANLVRSRQEGRWIYYSL---------"
      "-----------";
  seqs[counter++] =
      "---------------------------------"
      "QTTADIFKQLSDPTRIRIFWILCHCEECVINIASMMEMSSPAVAHHLRLLRSSGLIESRRDGKETYYRAVD-"
      "-----------------";
  seqs[counter++] =
      "-------------------------------------------"
      "ADEKRLKLVNLLLKQDYCVGALAKELEISKSAVSQHLKVLRESELVIGEKRGYWVHYSVQEDKLIEL-----"
      "-------";
  seqs[counter++] =
      "---------------------------------"
      "ERLTEIFKLLSDETRLRVVMLLAREETCVCEIVGVLGIPQPKVSKALSKLRDLGLVNDERKEKYVYY-----"
      "-----------------";
  seqs[counter++] =
      "------------------------------------"
      "AELFKALGHPLRLRILELLRTGEKTVGELQRLLMVEASSVSQQLAVMRAHHLVESRKQGTNVFYSVKD----"
      "--------------";
  seqs[counter++] =
      "------------------------------------"
      "AEFFRTLGHPARIRALELLSEREWSVSELVPEIGLEASHLSQQLGVLRRAGLVTTRKQGTTVFYAVASPEIV"
      "TL------------";
  seqs[counter++] =
      "------------------------"
      "MQTPATTIPHLIAAGFYALCDPLIISVLELLRQQELCVCDLCKALGVNQSKLSFHLKTLKETALVHSRQEGR"
      "WIY-----------------------";
  seqs[counter++] =
      "---------------------------------"
      "EELSQSFRVLGDPTRLRILRLVAEAPLNVTELVSLVGVAQSSVSHHLGKLKGLGLLREERHAGYSYYSL---"
      "-----------------";
  seqs[counter++] =
      "--------------------------------------------"
      "DEARLRLLVRLSEGERCVTDLAAGSDERMSTVSQRLKVLKGEGLVTGRREGKHVYYTLADRHV---------"
      "------";
  seqs[counter++] =
      "---------------------------MAVELVQSL----"
      "KALADDKRMQIIHLLLEGDLCVGALAQSLGISEPAVSQHLKVLREAGLVWGEKRG-----------------"
      "----------";
  seqs[counter++] =
      "-----------------------------------------------------"
      "MLFLKEYSVNEIAENLHLRQSTVSHQLRFLKNLRLVKYRREGTTLYYSHDDAHVMNMLKETINH-----";
  seqs[counter++] =
      "----------------------"
      "GAMAAAQPPIYRLKADFFRLLGHPARVRILELLRDGERAVGELQAALGLDSSGTSQHLTAMRRQGLLESRRA"
      "GTSVLYRVKDPRIFQLLEVA--------";
  seqs[counter++] =
      "------------------------------------"
      "AELLALLADRTRLALLHALTGGEADVSTLTQVCGAARPAVSQHLARLRLAGLVNTRKEGRRVIYSLRDGHLR"
      "RVVDEAL-------";
  seqs[counter++] =
      "-------------------------------------"
      "ETFRLIGDPSRLKILYILSHTEENVRNISAAFDMSPPAVSHHLRLLKSMKIIKSERRGKEVYYTL-------"
      "-------------";
  seqs[counter++] =
      "-------------------------------------"
      "KIYKVLSNMNRIKILYFLENHEADVSRIVDHVQLSQPIVSHQLAILYHYQLVTRHKRGKHVYYCLDDPHILE"
      "MVDAMLGHV----";
  seqs[counter++] =
      "---------------------------------"
      "KGVSQILKAIADENRAKIYALCQDEELCVCDIANILGVTIANASHHLRTLYKQGVVNFRKEGKLALYSLGDE"
      "HIRQIMMIALAHKKE--";
  seqs[counter++] =
      "------------------------------------"
      "ARVFKVLSVESRVRLIELLKQRSLCVNALARSLAITPAAVSQHLRVLRDAEVVIADKQGYHVHYRI------"
      "--------------";
  seqs[counter++] =
      "------------------------------"
      "KVAGQAAKLLAAIANARRLVILDIISQQETSVGSLAEQVGLSQSALSQHLAKLRSAKLVNTRRDAQTIYY--"
      "--------------------";
  seqs[counter++] =
      "-----------------------------------------"
      "VLANANRLLLMCQLSQGEKCVGELEELLDLHQPTLSQQLGVLRSEGLVSTRRDGKKIYYSVADARVLAL---"
      "---------";
  seqs[counter++] =
      "------------------------------"
      "EKSEQAARCLRAMAHPARLMILQLLSGSEMSVSELEKALDISQSNLSQHLNLMKDKQLLSSRRSGNQVYYSL"
      "KDPRLLGL------------";
  seqs[counter++] =
      "------------------------------------"
      "ANVFSLLSDPTRLRIILTLKEGEQPVGMIAEKLGRKPTIISQHLAKMRWGKLVRTRQEGTRIFYSLSDEHVS"
      "ALVDQAI-------";
  seqs[counter++] =
      "-------------------------------------"
      "DIFKALGDENRLRIINLLSKGKLCVCDIEAILMMTQSNVSRHLNKLKNVGIISSEKKSQWVYY---------"
      "-------------";
  seqs[counter++] =
      "-----------------------"
      "ELYQVEREELLSKAELLKVLGHPERLAIVLLTMDGERCVKELVEALGISQPKVSQHVGLMKELGILTFRKEG"
      "TKVLYRVNDRKVV--------------";
  seqs[counter++] =
      "-------------------------------------"
      "DIFRALGDPTRLRIVHLLRAMELAVGEIAQVVGQSQPRVSRHVRILAEAGLVERRKEGNWVFLRL-------"
      "-------------";
  seqs[counter++] =
      "---------------------------------------"
      "FHILQSDTRLRILFLLSQKQMCVCELEAGLDVTQSAISHSLSIMKNAGIVGVKREGRFAIYFIHDEEIRKMM"
      "QICRKYAEESR";
  seqs[counter++] =
      "---------------------------------------"
      "FGILSDKTRLRILLLLQNRELCVCEIFGALRMSQPRVSRQLAILKQSRIIKDRRSGKWIYYRIEE-------"
      "-----------";
  seqs[counter++] =
      "--------"
      "CPMIYALSRNIVISIFLMNISSAALQEIADFFEVLAVPTRLGILLAIGEREVCVCHLEAVLKLRQAAISQHL"
      "QVFKKNGWVISRRQGRFVYYKLSNPSVLPL------------";
  seqs[counter++] =
      "-----------------------------------"
      "MADIFDVVADPTRRDLLRVLPTGEISVSELVQTLGISQPTVSKHLRVLRDSGLVSVREEGQHRYYRLE----"
      "---------------";
  seqs[counter++] =
      "--------------------------------"
      "ALSATELFRLLGDETRLRAVVLLRRGELCVCELTETLGVSQPKMSRHLATLRDSGLVETRRSGQWIHYQL--"
      "------------------";
  seqs[counter++] =
      "----------------------------"
      "APDAPEQAAKFLKSLGHPDRLKVLCSLVGGEQSVASIEAQVGASQSAVSQHLSRLRSEGLLQARRDGRQVYY"
      "SIAD------------------";
  seqs[counter++] =
      "-----------------------------------------"
      "VLGNPDRLLLLCQLSQGEYAVGELETLLGITQPTLSQQLAVLREEQLVSTRREGKQVFYRIDSEAALALMQ-"
      "---------";
  seqs[counter++] =
      "---------------------------------------------------------------"
      "DLAQVLQMTPSAISHQLRVLKQMKLVTNRREGKTVFYSLADSHIKTIMNQGMEHIRK--";
  seqs[counter++] =
      "-------------------------------------"
      "EVLRVLADPTRLQLAGLLLDEEKSVSDLASQLDRPATGVSQHLAKMRMARLVSTRRRGTSVLYRVENDHVRQ"
      "LVVDTIGHVE---";
  seqs[counter++] =
      "------------------------------"
      "EVADRLAGIFKQVGDPTRLKIFWLLCQQEECVTNIAYLLDMSSPAISHHLKSLKLADLIESERKGKEMFY--"
      "--------------------";
  seqs[counter++] =
      "-------------------------------------"
      "DIFKQLSDPTRVRIFWLLSHREECVINIAALLDMSSPAVSHHLRSLTQSGLIESRRCGKEVYYKAGD-----"
      "-------------";
  seqs[counter++] =
      "-----------------------------------------"
      "VFAHPHRLMILSRLLRGECTVGEIDAATGIGQPALSQQLAQLRRAETVRTRREARQIHYSLADAHV------"
      "---------";
  seqs[counter++] =
      "-------------------------"
      "QLIGPELSRFKAEFFKALAHPLRIRIVDELRNGEVGVTHLCARLEVEQSSLSQQLAVLRARYIVNARKDGLS"
      "VLYSIRDPEIFSL------------";
  seqs[counter++] =
      "---------------------------------------"
      "FHALSDPIRLNILDILNNQEMCVGNICDLLSIKQSKVSFHLKILKESGFVETRQQGRCIYYRL---------"
      "-----------";
  seqs[counter++] =
      "------------------------------"
      "EVFESTARYFSVLGEPTRLKILHVICHKEKCVNDIIRATGLLQANVSRHLGLMYQAGLLSKRRDGTQIFYRV"
      "--------------------";
  seqs[counter++] =
      "-----------------------------------"
      "LADVFRLLGEPNRLRILCAIGSDCKSVSELMSETGIGQSNTSFHLRFLRNAALVNAEPRGRNMYYRVRDKEL"
      "LKL------------";
  seqs[counter++] =
      "--------------------------"
      "AAAKEELEEIASLLKLLGDKTRLTIFALLKVRELCVCELTELLHVSQPAISQHLRKLKLANLVRERKVGQWV"
      "HYSLRQRHIVLLEKSA--------";
  seqs[counter++] =
      "-------------------------------------"
      "ELFRILASQIKLEILSLLLENDLCVCQICAIVGTSQPNISQHLNTLRHLGVVDIRKDGTFIYYSL-------"
      "-------------";
  seqs[counter++] =
      "-----------------------"
      "EFERICPFMLETFETVAKAVADPSRVRILKLLEGGELCVCQITTVLDLAPATISKHLAALKTAGLVQQRRDG"
      "KWVYYRLAERDFNAYARSFLD------";
  seqs[counter++] =
      "--------------------------------------"
      "FLRAISDPNRLKILCVLQGGSKCVCEIVPLVGISDKLASHHLKQLKNVGLLIEKREGKFIRYNL-"
      "DKKVIKEYKNV--------";
  seqs[counter++] =
      "-----------------------------------"
      "LSDILHLMGEVNRLKLLIECLKGPKSVSDLAEQLQLSVPLTSHHLSLLRSARLLMANREGKHIYYSIYDAHV"
      "RCILEDMLKHFTE--";
  seqs[counter++] =
      "---------------------------------------"
      "FDLLSDPHRLELLSLHRAPGICVSDLAAALGRSENAVSQALRVLRQQGWVSSTRVGRAVSYRLDD-------"
      "-----------";
  seqs[counter++] =
      "--------------------------------"
      "AEEVSELLRILAHPERLMVLCQLTKGEVGVGQLQQSSALSQSAFSQHLTVLRKHGLIEARKESQQVFYSLAD"
      "TRVAQLIQ----------";
  seqs[counter++] =
      "--------------IEVEMATDEIMKKNAVEVA----"
      "ELLRVMAHPERLMVLCQLTHHEMGVGQLQQGSTLSQSAFSQHLTVLRKHGIIQARKESQQVFYRLADSRITA"
      "L------------";
  seqs[counter++] =
      "-----------------------------------------"
      "VLANPDRLKILCVLVDGEMNVQEIEESTDIHQPTLSQQLTVLRKADMVSTRREGKQIFYRLSDPKVLSLMQK"
      "LYEALNYC-";
  seqs[counter++] =
      "------------------------------------"
      "AEFFKTLGHPVRIRVLELLGQREHAVSEMLPEVGVEAANLSQQLAVLRRAGLVANRKEGSAVYYSL------"
      "--------------";
  seqs[counter++] =
      "------------------------------------"
      "ADLLLVMANAHRLRMLKTLAEREVAVNNLADIIGISQSALSQHLAKLRSRDLVKTRRDAQTIYY--------"
      "--------------";
  seqs[counter++] =
      "------------------------------"
      "EVFESVARYFSVLGEPTRLRILHALCQEEKCVNEIIKVTALAQANVSRHLGLMYQAGMLSRRREGTQIFYKV"
      "AD------------------";
  seqs[counter++] =
      "------------------------------"
      "ETFEKISDLFKQLGDPTRMRIFWILCHHEECVIHISARMDMSSPAVAHHLRLLKTSGLVTSRRQGKETYYRA"
      "SD------------------";
  seqs[counter++] =
      "------------------------------------------"
      "LSNEKRIRILYLLENHSFNVSELSEQLELPQPSVSHQLALLRQYQLVQAHRDGKQIFYTLDDPHIIEVLNDM"
      "LAHVQQ--";
  seqs[counter++] =
      "--------------------------------------"
      "FAKAISDPIRLRILYALREGELCVCELADALELRQSTLSTHLQIIRQAGLVQTRREGRWVYYALE-------"
      "------------";
  seqs[counter++] =
      "-------------------------------------"
      "EAFKAIADPTRRKILTLLRTGDLTAGEIASHFDMQKPSVSHHLKILKQADLVQDRREGQYIYYSL-------"
      "-------------";
  seqs[counter++] =
      "---------------------------------------"
      "FKALGQHLRLRIIALLAEQELCVCELEEILGITQPAISQHLRVLKEADLVWEEKVSQWVFYHLKKEKLAAVL"
      "QSWLAYLQ---";
  seqs[counter++] =
      "-----------------------------"
      "PKITGKWEDFFKVLSDETRLRILMLLNQRELCVCEICQILDLPQPKVSRHLAKMRDLDIVRGKKEDQWVFYY"
      "L--------------------";
  seqs[counter++] =
      "------------------------------"
      "ELFEEVANYFSLLCEPTRLKILYAVCNGERSVGDIVNEVESTQANVSRQINMLYRAKILARRKEGTQVYYRV"
      "DDEKTVDL------------";
  seqs[counter++] =
      "-------------------------------------"
      "EVFSMLADATRIRIILALRDQELSVNHLADIVDKSAPAVSQHLAKLRLARIVSTRQEGTKVFYRLTNEHARQ"
      "LVADAI-------";
  seqs[counter++] =
      "-------------------------------------"
      "ELLRALASPTRIAIVQSLGSESRCVHELVGELELSQPLVSQHLRVLKDAGVVRGERNGREIMYSLVDHHIVH"
      "IVDDALVHATE--";
  seqs[counter++] =
      "------------------------------"
      "ELLENAAATLRMLAEPTRLHLLWQLSQGPKSVTELTEAAAVPRTVVSQHLAKLRLSGMVDGRKNGRQVIYSL"
      "HDGHLVRLIRETINH-----";
  seqs[counter++] =
      "-------------------------------------"
      "EILKALSDENRLRILNLLRWGKLCVGEIQSILGITQSNASRHLNKLKGVGIIKFEKDAQWVHYKL-------"
      "-------------";
  seqs[counter++] =
      "------------------------------------------"
      "LSDPARLQMLWALSTEDLSLSDLAQLVGVSSTVASQLLSRLRTAGVLQTRKSGRHVIYSMHD----------"
      "--------";
  seqs[counter++] =
      "--------------------------------"
      "AEKTARMFKVLSVGSRVRMVELLKERSLCVNALARTLGITAAAVSQHLRVLRDAGLVCPEKHGYYVHYRI--"
      "------------------";
  seqs[counter++] =
      "---------------------------------"
      "ERLAEIFKALGHPTRVKIVEYLADGEKCVKDIWQEIGVPQPTVSQHINILKNAGIISFRKDG----------"
      "-----------------";
  seqs[counter++] =
      "------------------------"
      "MENLTPEAMEQVAAYFRALSEPTRLAILNLLREGERNVGELAQLCSCSPANVSRHLSLLSQHGLVRREGRGT"
      "AVYYRIADDSVYAL------------";
  seqs[counter++] =
      "------------------------------"
      "EKASDISKAFRHLGDPKRLQIFWLLCHRKECVINIAAIMGMSSPAISHHLKILKTAGLISSKREGKEMFYKA"
      "ND------------------";
  seqs[counter++] =
      "---------------------------------"
      "ENLSPLFHALADPNRLRIIELLRQEDLTVGSIAERLDISQPQTSKQLRVLYDAGLVS---------------"
      "-----------------";
  seqs[counter++] =
      "-------------------------------------"
      "ELLRALSAPIRLAIVSQLAEGERCVHELVNQLGAAQPLVSQHLRVLRGAGVVRGSRRGREIAYTLVDEHVAH"
      "IVADAVSHASE--";
  seqs[counter++] =
      "-------------------------------------"
      "ELLKIMAHPERMMVLCQLIEGEVAVAQLQQASLLSQSALSQHLALLRRQRLISARKRSQQVFYSLADQRVQQ"
      "LI-ASLQHIASC-";
  seqs[counter++] =
      "-----------------------------"
      "PEVFDRIAERLRILAHPHRLRMVEMLLAGKYSVGELAESCSIPSHMASEHLRLMQHCGLLGSEKEGRYTYYR"
      "I--------------------";
  seqs[counter++] =
      "------------------------------------------"
      "LSDQNRLRVLSLLDGNELTVKEMLEILQLSQSTLSSQLSQLKDSGLVQSRRDGQYVFYKLPRQYETQMVSNP"
      "ID------";
  seqs[counter++] =
      "------------------------------------"
      "AELLRQLANTNRLLILCHIAAEERSVGQLEADLGIKQPALSQQLAELRQYGLVKTRRQSRSIYYSIAD----"
      "--------------";
  seqs[counter++] =
      "------------------------------------------"
      "LSDETRMRILNLLEKGEMCVCEMEEILDISQSNASRHLTKLTNAEIINYNKVSKYVYYKI------------"
      "--------";
  seqs[counter++] =
      "--------------------------------------------------LLYQLKDGERCVGELV-"
      "VDGNKLSTVSARLQTLLNANLVKRRRDARHLYYRLADQHVVQLIDNALAHVDE--";
  seqs[counter++] =
      "--------------------------------"
      "AEEVAELLRVMAHPERLMVLCQLTQSEMGVGQLQQGSTLSQSAFSQHLTVLRKHGIIQARKESQQVFYRLAD"
      "SRI---------------";
  seqs[counter++] =
      "---------------------------------"
      "REMAELLGVLSHPCRVQIVEELRDSERNVNALQELLGISHSGVSQHLALLRTRKLLKERRSGRHVYYRL---"
      "-----------------";
  seqs[counter++] =
      "---------------------------------------"
      "FALLADPLRLRIVEALSREQLCTCHLVDITGARQTTISNHLRLLREAGVVASEPEGRYTWYRL---------"
      "-----------";
  seqs[counter++] =
      "-------------------------------------"
      "EILKALADETRIRILNLLYRETLCVCDLEEILKLSQSNASRHITKLKQAKLIAGEKQAQWIYYQV-------"
      "-------------";
  seqs[counter++] =
      "---------------------------------"
      "RNLVKFFAALADPTRLRLLNMMAGGEVCVCHFQGVLQTNQPKISRHLAYLKKAGLVEARRDGKWMHYRL---"
      "-----------------";
  seqs[counter++] =
      "---------------------------------------"
      "FKALGDEKRLRILSLLRQGERCACVLLEHLNLSQPTLSHHMKILCEARLVTGRKEGKWVYYSL---------"
      "-----------";
  seqs[counter++] =
      "-------------------------------------"
      "KLFKALAHPTRIQILNLLQEGELCVCEIYEALELSQSNISQHLKVLRDQNLVESQKVGVEVHYKIKN-----"
      "-------------";
  seqs[counter++] =
      "------------------------------------"
      "AEILSLLADRTRLALLRRLSLGEADVTTLTRACGVARPAVSQHLARLRLAGLVTTRKEGRRVVYALRHGHLR"
      "RLIDEALN------";
  seqs[counter++] =
      "-------------------------------------"
      "EIFKALSDKNRLLILDMISCGELCACDIMDVLNLTQPTISHHMKVLQKCELVDARKEGKWVFY---------"
      "-------------";
  seqs[counter++] =
      "-------------------------------------"
      "EVFRMLADATRVQVLWSLADREMSVNELAEQVGKPAPSVSQHLAKLRMARLVRTRRDGTTIFYRLENEHVRQ"
      "LVIDAVE------";
  seqs[counter++] =
      "---------------------------------------"
      "FKLLGNETRLNILLLLEKQPQTVSELVSALHLKQSNVSHQLAQLKHHQLIASTRRGKNLLYSLRDPHVITMI"
      "ETTYEH-----";
  seqs[counter++] =
      "---------------------------------------"
      "FKALADENRIRILNLLKNGKLCVCDIEAVLGIKQSNTSRHLNKLKMAGIIVSEKKSQWVYYRLND-------"
      "-----------";
  seqs[counter++] =
      "-------------------------------"
      "VYQLISEIFKTLAHPLRIQILMMLSEKERCVCELLNEIGVEQSNLSQHLRILKKQGIIDSRKDGQKMFYRI-"
      "-------------------";
  seqs[counter++] =
      "-----------------------"
      "ELPEMSSEQLARLASLFRLLGDEGRLKLVMACIDAPQPVCCLSEISGMSQPLTSHHLRGLREARILKSSRRG"
      "KQVLYELDDHHI---------------";
  seqs[counter++] =
      "------------------------------"
      "DMAQEQVTILKALADPNRLAIIQHLTEGEACVCELLQLFSVTQPTLSHHMRILSDADLVKGRREGKWIHY--"
      "--------------------";
  seqs[counter++] =
      "-------------------------------------"
      "DIFKALADENRIKIIKMLACCDMCVCDICGNLNLSQPAVSHHLKILSDSGLLNTTRKGKWIYYSL-------"
      "-------------";
  seqs[counter++] =
      "-------------------------------------"
      "DLFKALSDPTRRKILELLKEKDMSAGEIAEYFDISKPSISHHLNILKNAKLVLWEKDGQNIIY---------"
      "-------------";
  seqs[counter++] =
      "-------------------------------------"
      "KIFKILGSETRLNILLLLEKKDMTVTDLFNELEVSQPAISKQLAILKEYKIISYDKKGVENIYKLNDLHILN"
      "VINSTMGH-----";
  seqs[counter++] =
      "---------------------------"
      "MAEQVFAQVASYFGLLADPTRLRILSCLCGEERPVHDVVERIGLTQANISRHLNILYRAGVVDRRREGSSVL"
      "YKVVD------------------";
  seqs[counter++] =
      "------------------------------------"
      "AECLKALASPVRLKILFTLKDKPMCVTDLEQELGISQSSLSQHLRTLRYKGIVAKTRKGNKVYY--------"
      "--------------";
  seqs[counter++] =
      "------------------------------------"
      "ATLFKALAEPIRLRILALLKDGELCVCDLTETLALPQSTVSRHLAVLRTAGWIRGRKGGSWTYYSL------"
      "--------------";
  seqs[counter++] =
      "------------------------------------"
      "ADLFKALADPMRLRILALLRTREACVCELAGLLPITQPAVSQHLRKLRQAGLIHERRHKYWTYY--------"
      "--------------";
  seqs[counter++] =
      "------------------------------------"
      "ADFFKALAHPLRIRILEVLSEGERNVNELQTALGSEGSAVSQQLAVLRAKNLVNSFKEGTTVVYSLRD----"
      "--------------";
  seqs[counter++] =
      "-----------------------------"
      "PDRAGRIAEVLKAVAHPLRLRIVASLCREELNVSALAERLGASQAIVSQQLRILRSLGLVA-----------"
      "---------------------";
  seqs[counter++] =
      "-----------------------------------------------------"
      "VLKRRALCVTELTSQLGISQSATSQHLRVLKDARIVKFQKRGLHVYYHL--------------------";
  seqs[counter++] =
      "-----------------------------------"
      "LTNFLKIISDKNRLIILYLLSRNILCVCDIQKLIPLTQGALSIQLKNLMSAGLLESFKQGKWVFYKL-----"
      "---------------";
  seqs[counter++] =
      "-----------------------"
      "DMQMMMKDNANKASSLLKAISHESRLLILCLLLRREMTVGELAEYSSLSQSAFSQHLSVLRNNGLVKCRKEA"
      "QNVYYSIND------------------";
  seqs[counter++] =
      "--------------------------------------"
      "FFEAFSNKNRFEILMQLRNKELCAGELQQKLKIEQTNLSHDLKCLLNCRFISVRKDGR--------------"
      "------------";
  seqs[counter++] =
      "-----------------------------------------"
      "ILANEDRLLLLCQLSQGEKAVGELEDALGIHQPTLSQQLGVLRSDGLVNTRREGKRIFYSIADDKVLAL---"
      "---------";
  seqs[counter++] =
      "-------------------------------------------"
      "ADANRLRILACLKKGEVCVCDFTDFLNISQPAVSQHLRKLKEAGIITERKVGTWKHYRIQE-----------"
      "-------";
  seqs[counter++] =
      "-----------------------------------"
      "LSNIFKALNDPIRVKILFALLEYEICVGEMVNLLQIPQSHVSHQLRILRKYGIVEFTKDKKMSFYYIKNEYI"
      "KTL------------";
  seqs[counter++] =
      "---------------------------------------"
      "FDMLSAPNRLHLVWLLATGEFDVSTLAELSGSNVPAASQHLAKLRAAGIVTARRDGRRQLYRVEDPHIVTVI"
      "EQMFSHI----";
  seqs[counter++] =
      "--------------------"
      "KAKEAQLLSMEILEQAAECLRTLAHPHRLRIVQILLDHEESVGELARACELPSHMVSEHLRLLKDRGFLESR"
      "RDGRKVFY----------------------";
  seqs[counter++] =
      "-------------------------------------"
      "EIFKALGDENRIRILNLLIRQELCVCEIETVLDMTQSNASRHLNKLKTSGIITSEKKSQWVYYRV-------"
      "-------------";
  seqs[counter++] =
      "------------------------------"
      "EVEQYIDRFLDTVCDTRRRAIVELLAISEMRSGDIARAIGLSAATTSEHLRQLAQTGLLTSRRQGNTVYYSL"
      "CNHKLVQAFRDLLEAL----";
  seqs[counter++] =
      "------------------------------"
      "EMMEATARVLKLLGDPTRLTILAILQKRECCVCELMEVFSSSQPAISQHLRKLKDAGLLQEERRGQWVYYSL"
      "--------------------";
  seqs[counter++] =
      "------------------------------------"
      "AELLAVLGNERRLVILGHLTEGEISVGELAVLVGLSKSALSQHLSKLRKHQLVSTRRHRQTVYY--------"
      "--------------";
  seqs[counter++] =
      "-------------------------------------"
      "ELIRVLGDPLRLKIVTLLARETLCTSHLVEETGARQTNLSNHLRVLREAGVVETEPCGRFTYYKLRPDVIAA"
      "L------------";
  seqs[counter++] =
      "------------------------"
      "VSAFTEQFARKASDLLKAMSHETRLVILCLLSEKERSVGDIESILSMPQAAVSQQLARLRFDRLVKTRREGR"
      "TVYYSLASEEVTSL------------";
  seqs[counter++] =
      "------------------------------------"
      "ADLFKIMGDRSRLSMVAMMNRRECCVCDFTECFGMSQPAVSQHLKKLRAMGLIKERKEG-------------"
      "--------------";
  seqs[counter++] =
      "--------------------------------"
      "SQEAAKVMQLLSHPDRLLILCLLSEKEYSVGEIEKQLDIHQPMLSQHLNRLRQQSLVATRREGKYIYYQLCD"
      "------------------";
  seqs[counter++] =
      "--------------------------------"
      "ASAAAELLKLVANPNRLRILYLLTEGERSVSEIEQRLGIRQPTLSQQLGELRNAGTVTTRRAHKVVFYSL--"
      "------------------";
  seqs[counter++] =
      "--------------------------------------"
      "FFKALADDSRLKIVGILANQECSVEELAVLLQLKEPTVSHHLAKLKELNLVTMRPEGNSRLYQL--------"
      "------------";
  seqs[counter++] =
      "------------------------------------"
      "ASILKALGHPIRLKILYLLSEKEHCVCELLSQINTSQPNLSQHLSILRNLKLIKDERNGNMVIYKLQDNKIV"
      "--------------";
  seqs[counter++] =
      "-------------------------------------------------------------"
      "VGELTEEVGVSQSLVSQHLRLLRAGRLLKQTRSGRNVFYALPDCHVRTMLTNMMDHVLE--";
  seqs[counter++] =
      "------------------------------------------"
      "LASSNRLELLEALAQGERSVDALAQATGMSVANTSHHLQILRDSGLAESRKEGLQVIYRLSDDQIPVL----"
      "--------";
  seqs[counter++] =
      "---------------------------------"
      "QKLIKFFHALSDETRLKIIKLLEKSELCVCEIVAALDMVQPKVSFHLGVLKEAGLVKIKRKGKWILYSLDD-"
      "-----------------";
  seqs[counter++] =
      "--------------------------------"
      "AEVASELMKILSNENRLMILCQLVDGEKSVGELVELLDLNQPTVSQQLSRMKNQGLVSYRKNAQTVYYSL--"
      "------------------";
  seqs[counter++] =
      "-----------------------------------------"
      "LLGDKTRLTILSYLKDQELCVCELVDLLDISQPAISQHLKKLRVAEIIRERKQGTWVYYSL-----------"
      "---------";
  seqs[counter++] =
      "------------------------------------------"
      "LADSTRLKILNLLSRQEMAVCELIEALDLSQPAVSHHLKLLKQACLITDSREGKWVLY--------------"
      "--------";
  seqs[counter++] =
      "----------------------------------------------------------------"
      "MATEVGMEQSACSHQLRLLRNLGLVVGTRKGRSVVYSLYDNHVAELLDQAIYHIPVC-";
  seqs[counter++] =
      "--------------------------------------"
      "FLKVLGNPLRLQILKILSHVDMCVCAISEILGQQQTLVSHHLSKLKSARIVEERQNGKYRIYSIKDKRVKSL"
      "------------";
  seqs[counter++] =
      "--------------------------------"
      "AERLADRLRPLAQPQRLMILSLLLAGEHTVGEIETRTGIGQPALSQQLAELRRSGLVTTRRAARQVHYRIAD"
      "------------------";
  seqs[counter++] =
      "-------------------------------------"
      "DFLKLLADETRLKIIMMLSQRDMCVCEIMDELAMSQPAVSHHLRILKKSGIVRDDKDGRWVFYSL-------"
      "-------------";
  seqs[counter++] =
      "--------------------"
      "KTRELELSIPGVSDTLAKFFRAIADPNRLLLLEFLVSCEHTGNECVAHVRLAQSRVSSHLQCLVNCGFVRVR"
      "REGHFAYYRVVDERVIDL------------";
  seqs[counter++] =
      "------------------------------------------"
      "LADKNRLKIIQYLSTGQRNVSEVADRLNVEENLASHHLRVLASLGFLKNDKKGREVYYRINETRFVALLKDL"
      "L-------";
  seqs[counter++] =
      "---------------------SHELAAAAPGI-"
      "EAMAAVLALAGNEVRLKMLFLLLDQQLCVCDLADVLQMNVSAISQHLRKLKDGGVIQARKVGQTVFYSL---"
      "-----------------";
  seqs[counter++] =
      "------------------------------------------"
      "LADPTRMRMLWLISGEEYDVASLAAAVDIARPAVSQHLAKLKLAGLVTQRRDGRRILYRARGGHVLAEVMNA"
      "ADH-----";
  seqs[counter++] =
      "--------------------------------"
      "AREVSRLLSVLANENRLLIVCLMMRSEMKVGELVDALHLSQSALSQHLTKLREEGLVEFRRESQTLHYKIAD"
      "ERVTKL------------";
  seqs[counter++] =
      "-----------------------------------"
      "LTNIFKVLSDENRLRMIVLLYQEELCVCELSGILNVPQPRISQNLSRLRDLNLVDDERKEKFVFYSL-----"
      "---------------";
  seqs[counter++] =
      "------------------------------------"
      "AEILRILSHPERLLVLCQLMEGELGAGQLQNSSTLSQSAFSQHLTVLRKHNLVKVRKESQQVFYSLADERIA"
      "ALIHN---------";
  seqs[counter++] =
      "------------------------------------------"
      "MGNPQRLRILLLLAEHERSVIELEALVGLSQSAVSQHLARLRQIKLVRFRRDGQMTFYAL------------"
      "--------";
  seqs[counter++] =
      "-----------------------------PIYAQ-"
      "LARVGKALASPIRLRLLDLLDGAELTVEELSEQAGVPLKNTSAQLQQLRAANLVATRKEGTRVHYRLAD---"
      "---------------";
  seqs[counter++] =
      "---------------------------"
      "MATDALDQVSHLFKLMGHPKRLQLLYLLIQQSMTVSQISERLKWEQSAVSHQLQVLRKYQIVERVKNGRQVV"
      "YRLVD------------------";
  seqs[counter++] =
      "------------------------------------------"
      "LADSTRVQVLWALVDRELSVNDLAEHVGKPAPSVSQHLAKLRMARLVRTRKEGTQVIYRLENDHVRQLVTDA"
      "VN------";
  seqs[counter++] =
      "------------------------------------"
      "AEFFKTLGHPVRIRVLELLSEREHAVSEMLNEVGVEAAHLSQQLAVLRRAGLVTARREGSAVHYTLAD----"
      "--------------";
  seqs[counter++] =
      "------------------------------------------"
      "LADENRLRILRALVGTEKPVSKLVEELGISQPLVSHHLKELRRALLVSVERRGPFVYCRLAD----------"
      "--------";
  seqs[counter++] =
      "--------------------------------"
      "AEHVAEMLKLMAHPHRLMILCLLVESEHNVGELVEALDINQTALSNHLSKLRSAGLIDYTRYHRVLQYRL--"
      "------------------";
  seqs[counter++] =
      "---------------------------"
      "IVPVVFQGAADLFAALSCPTRLRIVCALSQADHTVRDLARASQCSQANVSGHLRLLRRANIVRCERSGNYVL"
      "YHL--------------------";
  seqs[counter++] =
      "------------------------------"
      "EVFEEVANYFSLLCEPTRLKILYAVCNGERSVGDIVAQVESTQANVSRQLAMLYRAKILARRKEGTLVFYRV"
      "DD------------------";
  seqs[counter++] =
      "---------------------------------"
      "RQLADVGGALSNPHRLKMISLLAQGDKPIDELAKLTNQSLAAASANVKVLRNCHLIATEKRGRSVYCSLKDP"
      "RVAELW-----------";
  seqs[counter++] =
      "-------------------------------------"
      "ELLRALSAPIRLAIVSELAEGERCVHELVDKLGAPQPLVSQHLRVLRSAGVVRGSRRGREIAYTLVDEHVAH"
      "IVTDAVSH-----";
  seqs[counter++] =
      "-------------------------------------------"
      "AHPLRLKILCVLGEGEACVQDIVEAVGTSQSNISQHLAILRDKGVLQTRKDANRVYYRVGDQRTLQL-----"
      "-------";
  seqs[counter++] =
      "-----------------------------------------"
      "VAAEPTRRRLLQLLAPGERTVTQLASQFTVTRSAISQHLGMLAEAGLVTARKQGRERYYRL-----------"
      "---------";
  seqs[counter++] =
      "-----------------------------------------"
      "VAAEPTRRRLLQLLAPGERTVTQLASQFTVTRSAISQHLGMLAEAGLVTARKQGRERYYRL-----------"
      "---------";
  seqs[counter++] =
      "-------------------------------------"
      "EMLKILSDTNRLRILNLLYIQELCVCELEYLLTISQSNLSKHLRLMGEIGFLDSRRQNKFIYYKI-------"
      "-------------";
  seqs[counter++] =
      "------------------------------------"
      "AESFRLLADPTRIKVLWALLQGESSVACLAELAGAAPTXXSXHXXKLRLAGLVTGRREGTFVYYSAVNNHVR"
      "GLLAQALFH-----";
  seqs[counter++] =
      "--------------------------------AEATATLFA-"
      "LANQNRLLLLCQLCNGEMSVSALEEALGIHQPTLSQQLGVLRSEGLIASRREGKRIYYSVANPKVLVLINTL"
      "VD------";
  seqs[counter++] =
      "------------------------------"
      "EILDAAGELLRALAAPVRIAIVLQLRESQRCVHELVDALHVPQPLVSQHLKILKAAGVVTGERSGREVLYRL"
      "ADHHLAHIVLDAVAH-----";
  seqs[counter++] =
      "---------------------------------"
      "QDLLKVFYALSDSVRLGIVSLLECEELCVCQITQAFGLSQPNASFHLRVLREANLVLWEKRGKWTYYKINHH"
      "-----------------";
  seqs[counter++] =
      "---------------------------------------"
      "FKALGEPSRLKIIKLLSQQSMCVCELSEVLDMSQPRVSQHLRTLKEVDLVYEERQGFWTYYKL---------"
      "-----------";
  seqs[counter++] =
      "---------------------------------------"
      "FKALSDQTRLRMVTLLSRREYCNCEFVSIFGISQPAISRHIARLKEARLIHERRPGQWIYYSL---------"
      "-----------";
  seqs[counter++] =
      "---------------------------------------"
      "FKALADSNRLRILDYLKKGKSCACDLSDNLGIPQTALSYHMRILCQAKLVKSEQVGKWKHYQLND-------"
      "-----------";
  seqs[counter++] =
      "-----------------------ELEAKAEDAAQ----"
      "FLKMIASPPRLLLLCHMAERECSVGELAERTGMRMPTVSQQLSLLRAQGLVNTRRDGTTIYYRL--------"
      "------------";
  seqs[counter++] =
      "-------------------------------------"
      "KIFSALSDKNRLRIYLLLTQAELCVCELVNILDMEQSRISHSLRILKEAKLINNHRVGK-------------"
      "-------------";
  seqs[counter++] =
      "-----------------------------"
      "PESLRSIASLLKALSDPLRLQVLEQLSTGERCVCDLTSSLALSQSRLSFHLKVMKEAGLLSDRQSGRWVYYR"
      "IRPESLNAL------------";
  seqs[counter++] =
      "------------------------------------"
      "AKLFRGFADPSRLAILEVLRSAPATVGEIAASTGQGLSNVSNHLRCLRDCGLVVRQRDGQRVRYSLSDQRVA"
      "AL------------";
  seqs[counter++] =
      "--------------------------AMAGNVEQA-"
      "EQLLKVLANKNRLMILCSLQDSEMSVSQLNEAVPLAQSALSQHLAALRKANIVATRRESQTIYYRVIDENAV"
      "VL------------";
  seqs[counter++] =
      "-------------------------------------------"
      "ADTLRVQVLSLLHKNSFSVGELVEILGVRQSALSHHLKVLAQAELVATRREGNSIFY---------------"
      "-------";
  seqs[counter++] =
      "----------------------------"
      "AAQVFERAAELFGLLSSPLRLRIVGELCRGELNVGQLQERIGATQSNMSQHLSVLYRAGVVARRRDGAQVHY"
      "RI--------------------";
  seqs[counter++] =
      "----------------------------------"
      "ALLRWLKALSDDTRLRLLHLLSRYELSVGEVVQVLGMSQPRVSRHLKILADAGMVQVRRDGLWAFYSATSH-"
      "----------------";
  seqs[counter++] =
      "-----------------------"
      "DMQDIAAQLQELHARVCKAIADPKRLLIINELRDRELSVGELCEATGLSQSNASQHLTILRERGIVTTRRVK"
      "NNVFYSLRSQKIV----QAVDLLRE--";
  seqs[counter++] =
      "---------------------------------------"
      "FRALSDRTRLRILNLLRGGELCVCDLVDVLDVPQPTASRHLAYLRNAGLVLARKEGLWHYYRL---------"
      "-----------";
  seqs[counter++] =
      "------------------------------------"
      "AAYFQALAEPTRLQILNFLRQQERNVGELAQLCGYSSANISRHLALLTQHGLVSRQARGNSAYYRIAD----"
      "--------------";
  seqs[counter++] =
      "-------------------------------------"
      "DMFKAMADPTRRRILQLLSEKNLSAGEIAEEFTMSKPAISKHLDILKTSELITCEKQGQYVIYAINTSAVEQ"
      "MYCRFLD------";
  seqs[counter++] =
      "--------------------------------"
      "AHEASDLLKALAHQTRLLILCILANEERTVGEIENILGIQQAMVSQQLARLRLEGLVHTRRQGRLVYYSIGN"
      "VSVLAFLESLFD------";
  seqs[counter++] =
      "--------------------------------------------"
      "DPLRLNVLRALANDSFGVLELAQIFGIGQSGMSHHLKVLAQADLVATRREGNAIFYRRALPHLNAVKDGALD"
      "H-----";
  seqs[counter++] =
      "------------------------------"
      "EMTPDIVSFLKTISEENRLKILCFLRDWEKCVCEIVEFLKIPQNLVSHHLRKLKDARILSARKDGMNVRYSI"
      "NEDEI---------------";
  seqs[counter++] =
      "----------------------------"
      "AEPVFDRLAQVLDLAGNANRLKIIYLLEESNLCVCDLSDILGMSIPAVSQHLRKLKDAQLIQARKVGQTVFY"
      "SL--------------------";
  seqs[counter++] =
      "--------------------------------"
      "ASAAARLMKLMANEQRLILMCRLGEGECSVGDLAAHVGLAQSAASQHLAKLRAEGVVATRRDGQTIYYRLED"
      "------------------";
  seqs[counter++] =
      "------------------------------"
      "ENATEVAGILKQLSNPYRLMILCCLSENELTVGDLNQRIDLSQSALSQHLAKLRESNIVTTRRESQTIFYRI"
      "--------------------";
  seqs[counter++] =
      "-----------------------------------------"
      "VLSNPDRLKILCVLIDGELNVQQIEKTAQVYQPTLSQQLTILRKSKMVSTRREGKQIFYQFSDMRILQIMQT"
      "LYD------";
  seqs[counter++] =
      "---------------------------------------"
      "FHSLSDATRLAIVLRLARGEARVADLVGELGLAQSTVSAHVACLRDCQLVAGRPEGRQIFYRLARRELIDLL"
      "ASALE------";
  seqs[counter++] =
      "-----------------------------------"
      "MAKFYRALGDPTRLDLLEFCAEDERTGNECVERAGLSQGRVSAHLACLVSCGLVSVRRQGRFAYYRVTDPRV"
      "AEL------------";
  seqs[counter++] =
      "---------------------------"
      "IETEISPELTNFIKVLSNPIRAGIIKMLKKRWMCVCLIAKALNQDQTLISHHLRTLKNMNLLHERREGK---"
      "-----------------------";
  seqs[counter++] =
      "------------------------------------"
      "AQVFKALGNPVRMALVQELLAGERCVCDLAQALGGNMPAVSKHLATLREAGIVSCRREGTTIHYSL------"
      "--------------";
  seqs[counter++] =
      "------------------------------------"
      "AKLMEMLSQPVRLRILCILLDGEQSVLKLADMAGLSQPAMSHHLRKLRDADLVNTRRDAQTIYYSLKGQEVS"
      "AV------------";
  seqs[counter++] =
      "------------------------------------------"
      "LGSSNRLMLVCQLLDGERAVGELAEALGLAQSVVSQHLSLLRRDGLVTGRRDGQSIYYAISDDRVHAL----"
      "--------";
  seqs[counter++] =
      "----------------------------------------"
      "AVLANINRLLLMCQLSQGEKCVGELEELLDLHQPTLSQQLGVLRGAGLVNTRRDGKKIHYSVADARVLTL--"
      "----------";
  seqs[counter++] =
      "---------------------------------"
      "QTLLGFFQALADANRLRIVGVLAQGPQTVEQISALLGLGMSTTSHHLRKLAKAGLVEARADGHYSVYSLRTQ"
      "TLEELAKNLL-------";
  seqs[counter++] =
      "-------------------------------------"
      "DLFKCIGNPTRYKILKVLCERPLCVNKLNEAVGYSQPNISQHLKLMRMSGIVTCSKNGMNICYQIADDDIIK"
      "LLELAEDILKNRR";

  char** seqsCpy = new char*[counter];
  for (int k = 0; k < counter; ++k) {
    seqsCpy[k] = MultipleAlignment::initX(122);
    for (int pos = 0; pos < 122; ++pos) {
      //            seqs[k][pos] = (seqs[k][pos] == '-') ?
      //            MultipleAlignment::GAP : subMat.aa2num[(int) seqs[k][pos]];
      seqsCpy[k][pos] =
          (seqs[k][pos] == '-')
              ? MultipleAlignment::GAP
              : static_cast<int>(subMat.aa2num[(int)seqs[k][pos]]);
    }
  }

  MultipleAlignment::MSAResult res(122, 122, counter, seqsCpy);
  MultipleAlignment::print(res, &subMat);

  MsaFilter msaFilter(10000, counter, &subMat, par.gapOpen.aminoacids,
                      par.gapExtend.aminoacids);
  std::vector<Matcher::result_t> empty;
  size_t filteredSetSize = msaFilter.filter(res, empty, 0, 0, -20.0, 90, 100);

  /*    std::cout << "Filtered:" << filterResult.setSize << std::endl;
  //    for(size_t k = 0; k < res.setSize; k++){
  //        std::cout << "k=" << k << "\t" << (int)filterResult.keep[k] <<
  std::endl;
  //    }

      std::cout <<"Filterted MSA" << std::endl;
      for(size_t k = 0; k < filterResult.setSize; k++){
          printf("k=%.3zu ", k);
          for(size_t pos = 0; pos < res.centerLength; pos++){
              char aa = filterResult.filteredMsaSequence[k][pos];
              printf("%c", (aa < MultipleAlignment::NAA) ?
  subMat.num2aa[(int)aa] : '-' );
          }
          printf("\n");
      }
  */
  // seqSet.push_back(s5);
  PSSMCalculator pssm(&subMat, par.maxSeqLen, filteredSetSize, 1.0, 1.5);
  PSSMCalculator::Profile profile = pssm.computePSSMFromMSA(
      filteredSetSize, res.centerLength, (const char**)res.msaSequence, false);
  // std::string libraryString((const char *)Library1_lib, Library1_lib_len);
  ProfileStates ps(32, subMat.pBack);

  char mmOrder[] = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M',
                    'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y', '\0'};
  for (size_t i = 0; i < 20; i++)
    std::cout << mmOrder[ProfileStates::hh2mmseqsAAorder(i)] << std::endl;

  //    std::ifstream libFile;
  //    libFile.open ("/home/clovis/Software/mmseqs-dev/data/LibraryPure.lib");

  //    std::string libData((std::istreambuf_iterator<char>(libFile)),
  //                        std::istreambuf_iterator<char>());

  //    ps.read(libData);

  //    std::string
  //    sequence({0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,'\0'});

  size_t L = res.centerLength;
  std::string disc;
  ps.discretize(profile.prob, L, disc);
  for (size_t i = 0; i < L; i++)
    std::cout << i << "," << (int)disc[i] << std::endl;
  // pssm.printProfile(res.centerLength);
  pssm.printPSSM(res.centerLength);
  for (int k = 0; k < counter; ++k) {
    free(seqsCpy[k]);
  }
  delete[] seqsCpy;
  return 0;
}

// PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNF
//                     ALLDTGADDTVISEEDWPTDWPVMEAANPQIHGIGGGIPVRKSRDMIELGVINRDGSLERPLLLFPLVAMTPVNILGRDCLQGLGLRLTNL
