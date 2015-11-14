<!--
************************************************************************

              WebTrilinos: A Web Interface to Trilinos
                 Copyright (2006) Sandia Corporation

Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
license for use of this work by or on behalf of the U.S. Government.

This library is free software; you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as
published by the Free Software Foundation; either version 2.1 of the
License, or (at your option) any later version.
 
This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.
 
You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301
USA

Questions about WebTrilinos? Contact Marzio Sala (marzio.sala _AT_ gmail.com)

************************************************************************
-->

<?php 

  $styleSheetString = "";
  $altStyleSheetString = "";
  $filename = "";
  $title = "";
  $styles = "";
  $scripts = "";
  $bodyAttributes = "";
  $dir = "";

  function setFilename ($in_filename) { global $filename; $filename = $in_filename; }   
  function setTitle ($in_title) { global $title; $title = $in_title; }    
  function setStyles ($in_styles) { global $styles; $styles = $in_styles; }    
  function setScripts ($in_scripts) { global $scripts; $scripts = $in_scripts; }  
  function setBodyAttributes ($in_bodyAttributes) { global $bodyAttributes; $bodyAttributes = $in_bodyAttributes; }  
  function setDir ($in_dir) { global $dir; $dir = $in_dir; } 
  
  function includeStyleSheet ($in_styleSheet) { 
    global $styleSheetString; 
    $styleSheetString .= "<link rel=\"stylesheet\" type=\"text/css\" href=\"$in_styleSheet\" \\>\n"; 
  }
  
  function includeAltStyleSheet ($in_altStyleSheet, $in_title) { 
    global $altStyleSheetString; 
    $altStyleSheetString .= "<link rel=\"alternate stylesheet\" type=\"text/css\" href=\"$in_altStyleSheet\" title=\"$in_title\" \\>\n"; 
  }
  
?>
