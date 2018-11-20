<!DOCTYPE html>
<html xmlns="http://www.w3.org/1999/xhtml">
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1.0" />
    <title>Free Bootstrap Admin Template : Two Page</title>
    <!-- BOOTSTRAP STYLES-->
    <link href="../assets/css/bootstrap.css" rel="stylesheet" />
    <!-- FONTAWESOME STYLES-->
    <link href="../assets/css/font-awesome.css" rel="stylesheet" />
    <!-- CUSTOM STYLES-->
    <link href="../assets/css/custom.css" rel="stylesheet" />
    <link href="../assets/css/style.css" rel="stylesheet" />
  </head>
  <body>
    <div class="navbar" style="background:#DFDFDF; border: 2px solid #DFDFDF; font-size:120%;">
      <ul class="nav navbar-nav navbar-left" style="height=120px;">
	<li><a href="../index.html"><font color="black" size=4px>TCGAportal</font></a></li>
	<li><a href="../index.html">TCGA</a></li>
	<li><a href="index.html">Pan-cancer</a></li>
	<li><a href="../download.html">Download</a></li>
	<li><a href="../news.html">News</a></li>
	<li><a href="../about.html">About</a></li>
      </ul> 
    </div>

    <?php
      //<!--JS 页面自动刷新 auto refresh page-->
      echo ("<script type=\"text/javascript\">");
      echo ("function fresh_page()");    
      echo ("{");
      echo ("window.location.reload();");
      echo ("}"); 
      echo ("setTimeout('fre搜索sh_page()',1000);");      
      echo ("</script>");

      //-----------------------------------------------------------
      $file = '../../visitor.log';
      $current = file_get_contents($file);
      $mydate=getdate(date("U"));
      $current .= $_SERVER['REMOTE_ADDR']."\t"."$mydate[weekday], $mydate[month] $mydate[mday], $mydate[year].\t"."Pan-TCGA"."\t".$_GET["geneName"]."\n";
      file_put_contents($file, $current);

      $GeneName=$_GET["geneName"];
      echo "<strong>";
      echo "<p></p>";
      echo "<p>";
      echo $GeneName;
      echo "</p>";
      echo "</strong>";
      echo '<hr/>';
      $figure=glob("TCGA/figures/".$GeneName."_*.png");
      if( $GeneName!="" && !empty($figure) ){

      echo 'Expression in TCGA','<br/>';
      echo '<img src='.$figure[0].' width=100% heigth=100%>';
      echo '<br/>';
      echo '<hr/>';
      
      echo 'Expression in GTEx normal tissue','<br/>';
      $figure=glob("GTEx/all_gene_distribution/".$GeneName."_*log2.png");
      echo '<img src='.$figure[0].' width=100% heigth=100%>';
      echo '<hr/>';
      
      }else{
      echo "can not find ";
      echo "<strong>";
      echo $GeneName." ";
      echo "</strong>";
      }
      //-----------------------------------------------------------

      ?>

  </body>
</html>
