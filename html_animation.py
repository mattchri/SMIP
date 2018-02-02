# Description: Script used to generate an animated html and gif from a set of
#images. The code creates an html script in the directory where the png files
#are located. The html code is constructed from an html_template as described 
#below. This code simply inserts new lines of text containing the 
#pre-generated png files into the new html file.
#
# Required Inputs:
#imagefiles: array of strings containing each png filename
#variablename: name of variable at START of each filename
#outPath: output path of script
#
# OPTIONAL Input:
#make_gif: animate images into a gif
#
#02/01/18, MC: upload initial version of html_animation.py to the repo
#------------------------------------------------------------------------------

#import libs
import os

# Top-level code: html_script
def html_script( imagefiles, variablename, outPath, make_gif = None):
    #read html template
    lines = html_template()
    lct = len(lines)

    #old way by reading ascii template file (see notes below)
    #lines = [line.strip() for line in open(template, 'r')]

    #number of imagefiles to import to html file
    fct = len(imagefiles)

    print('hi')
    print(make_gif)
    print('hi')

    #format new lines to add to HTML file script
    newLines = []
    for i in range(fct):
        if i+1 < 10:
            istr=str(i).zfill(1)
        if i+1 >= 10 and i+1 < 100:
            istr=str(i).zfill(2)
        if i+1 >= 100 and i+1 < 1000:
            istr=str(i).zfill(3)
        newLines.append('modImages['+istr+'] = '+'"'+'./'+imagefiles[i]+'";')
    if fct < 10:
        nstr = str(fct).zfill(1)
    if fct >= 10 and fct < 100:
        nstr = str(fct).zfill(2)
    if fct >= 100 and fct < 1000.:
        nstr = str(fct).zfill(3)
    
    #create html script: add new lines to template
    fname=outPath+'animate_'+variablename+'.html'
    lines_for_file = []
    print('creating: '+fname)
    fout1 = open(fname,"w")
    for k in range(lct):
        if lines[k] != 'last_image = 96;':
            fout1.writelines( lines[k]+"\n")
            lines_for_file.append(lines[k])

        if lines[k] == 'last_image = 96;':
            fout1.writelines('last_image = '+nstr+';' + "\n")
            lines_for_file.append('last_image = '+nstr+';')

        if lines[k] == 'modImages = new Array();':
            fout1.writelines( '' + "\n")
            lines_for_file.append('')
            for ik in range(fct):
                fout1.writelines( newLines[ik] + "\n")
                lines_for_file.append(newLines[ik])
    fout1.close()

    #Another way to write a list of srtings to a file
    #fout1 = open(fname,"w")
    #fout1.writelines( "%s\n" % item for item in lines_for_file)
    #fout1.close()

    #create gif
    if make_gif == 'yes':
        wildTSTR = outPath + variablename + '*' + '*.png'
        gifFile  = outPath + 'animate_' + variablename + '.gif'
        print('creating: ',gifFile)
        os.system('convert -delay 15 -loop 0 '+wildTSTR+' '+gifFile)

    print('end of html_script')
    errval=0
    return errval

#---------------------------------------------------------------------------
# Template HTML file was converted to a list
#Original file can be found here: 
#http://storm.colorado.edu/~cassano/atoc4720/Daily_images/20171025/gom_vis.html
#Or on CEMS: '/home/users/mchristensen/python/animate_template.html.orig'
#---------------------------------------------------------------------------
# Function containing the HTML template - output has been converted to a python list
def html_template():
    template_lines = ['<HTML>', '', '<HEAD>', '', '<TITLE> Animation using Javascript Image Player </TITLE>', '', '<SCRIPT LANGUAGE="JavaScript">', '// <!--', '', '//============================================================', '//                >> jsImagePlayer 1.0 <<', '//            for Netscape3.0+, September 1996', '//============================================================', '//                  by (c)BASTaRT 1996', '//             Praha, Czech Republic, Europe', '//', '// feel free to copy and use as long as the credits are given', '//          by having this header in the code', '//', '//          contact: xholecko@sgi.felk.cvut.cz', '//          http://sgi.felk.cvut.cz/~xholecko', '//', '//============================================================', '// Thanx to Karel & Martin for beta testing and suggestions!', '//============================================================', '//', '//     modified by D. Watson and A. Earnhart (CIRA/CSU), 7/30/97', '//     and Greg Thompson (NCAR/RAP) May 07, 2000', '//', '//============================================================', '', '//********* SET UP THESE VARIABLES - MUST BE CORRECT!!!*********************', '', 'modImages = new Array();', '', '', '', 'first_image = 1;', 'last_image = 96;', '', '//**************************************************************************', '', '//=== THE CODE STARTS HERE - no need to change anything below ===', '', '//=== global variables ====', 'theImages = new Array();      //holds the images', 'imageNum = new Array();       //keeps track of which images to omit from loop', 'normal_delay = 150;', 'delay = normal_delay;         //delay between frames in 1/100 seconds', 'delay_step = 50;', 'delay_max = 6000;', 'delay_min = 10;', 'dwell_multipler = 5;', 'dwell_step = 1;', 'end_dwell_multipler   = dwell_multipler;', 'start_dwell_multipler = dwell_multipler;', 'current_image = first_image;     //number of the current image', 'timeID = null;', 'status = 0;                      // 0-stopped, 1-playing', 'play_mode = 0;                   // 0-normal, 1-loop, 2-sweep', 'size_valid = 0;', '', '//===> Make sure the first image number is not bigger than the last image number', 'if (first_image > last_image)', '{', 'var help = last_image;', 'last_image = first_image;', 'first_image = help;', '}', '', '//===> Preload the first image (while page is downloading)', 'theImages[0] = new Image();', 'theImages[0].src = modImages[0];', 'imageNum[0] = true;', '', '//==============================================================', '//== All previous statements are performed as the page loads. ==', '//== The following functions are also defined at this time.   ==', '//==============================================================', '', '//===> Stop the animation', 'function stop()', '{', '//== cancel animation (timeID holds the expression which calls the fwd or bkwd function) ==', 'if (status == 1)', 'clearTimeout (timeID);', 'status = 0;', '}', '', '', '//===> Display animation in fwd direction in either loop or sweep mode', 'function animate_fwd()', '{', 'current_image++;                      //increment image number', '', '//== check if current image has exceeded loop bound ==', 'if (current_image > last_image) {', 'if (play_mode == 1) {              //fwd loop mode - skip to first image', 'current_image = first_image;', '}', 'if (play_mode == 2) {              //sweep mode - change directions (go bkwd)', 'current_image = last_image;', 'animate_rev();', 'return;', '}', '}', '', '//== check to ensure that current image has not been deselected from the loop ==', "//== if it has, then find the next image that hasn't been ==", 'while (imageNum[current_image-first_image] == false) {', 'current_image++;', 'if (current_image > last_image) {', 'if (play_mode == 1)', 'current_image = first_image;', 'if (play_mode == 2) {', 'current_image = last_image;', 'animate_rev();', 'return;', '}', '}', '}', '', 'document.animation.src = theImages[current_image-first_image].src;   //display image onto screen', 'document.control_form.frame_nr.value = current_image;                //display image number', '', 'delay_time = delay;', 'if ( current_image == first_image) delay_time = start_dwell_multipler*delay;', 'if (current_image == last_image)   delay_time = end_dwell_multipler*delay;', '', '//== call "animate_fwd()" again after a set time (delay_time) has elapsed ==', 'timeID = setTimeout("animate_fwd()", delay_time);', '}', '', '', '//===> Display animation in reverse direction', 'function animate_rev()', '{', 'current_image--;                      //decrement image number', '', '//== check if image number is before lower loop bound ==', 'if (current_image < first_image) {', 'if (play_mode == 1) {               //rev loop mode - skip to last image', 'current_image = last_image;', '}', 'if (play_mode == 2) {', 'current_image = first_image;     //sweep mode - change directions (go fwd)', 'animate_fwd();', 'return;', '}', '}', '', '//== check to ensure that current image has not been deselected from the loop ==', "//== if it has, then find the next image that hasn't been ==", 'while (imageNum[current_image-first_image] == false) {', 'current_image--;', 'if (current_image < first_image) {', 'if (play_mode == 1)', 'current_image = last_image;', 'if (play_mode == 2) {', 'current_image = first_image;', 'animate_fwd();', 'return;', '}', '}', '}', '', 'document.animation.src = theImages[current_image-first_image].src;   //display image onto screen', 'document.control_form.frame_nr.value = current_image;                //display image number', '', 'delay_time = delay;', 'if ( current_image == first_image) delay_time = start_dwell_multipler*delay;', 'if (current_image == last_image)   delay_time = end_dwell_multipler*delay;', '', '//== call "animate_rev()" again after a set amount of time (delay_time) has elapsed ==', 'timeID = setTimeout("animate_rev()", delay_time);', '}', '', '', '//===> Changes playing speed by adding to or substracting from the delay between frames', 'function change_speed(dv)', '{', 'delay+=dv;', '//== check to ensure max and min delay constraints have not been crossed ==', 'if(delay > delay_max) delay = delay_max;', 'if(delay < delay_min) delay = delay_min;', '}', '', '//===> functions that changed the dwell rates.', 'function change_end_dwell(dv) {', 'end_dwell_multipler+=dv;', 'if ( end_dwell_multipler < 1 ) end_dwell_multipler = 0;', '}', '', 'function change_start_dwell(dv) {', 'start_dwell_multipler+=dv;', 'if ( start_dwell_multipler < 1 ) start_dwell_multipler = 0;', '}', '', '//===> Increment to next image', 'function incrementImage(number)', '{', 'stop();', '', '//== if image is last in loop, increment to first image ==', 'if (number > last_image) number = first_image;', '', '//== check to ensure that image has not been deselected from loop ==', 'while (imageNum[number-first_image] == false) {', 'number++;', 'if (number > last_image) number = first_image;', '}', '', 'current_image = number;', 'document.animation.src = theImages[current_image-first_image].src;   //display image', 'document.control_form.frame_nr.value = current_image;                //display image number', '}', '', '//===> Decrement to next image', 'function decrementImage(number)', '{', 'stop();', '', '//== if image is first in loop, decrement to last image ==', 'if (number < first_image) number = last_image;', '', '//== check to ensure that image has not been deselected from loop ==', 'while (imageNum[number-first_image] == false) {', 'number--;', 'if (number < first_image) number = last_image;', '}', '', 'current_image = number;', 'document.animation.src = theImages[current_image-first_image].src;   //display image', 'document.control_form.frame_nr.value = current_image;                //display image number', '}', '', '//===> "Play forward"', 'function fwd()', '{', 'stop();', 'status = 1;', 'play_mode = 1;', 'animate_fwd();', '}', '', '//===> "Play reverse"', 'function rev()', '{', 'stop();', 'status = 1;', 'play_mode = 1;', 'animate_rev();', '}', '', '//===> "play sweep"', 'function sweep() {', 'stop();', 'status = 1;', 'play_mode = 2;', 'animate_fwd();', '}', '', '//===> Change play mode (normal, loop, swing)', 'function change_mode(mode)', '{', 'play_mode = mode;', '}', '', "//===> Load and initialize everything once page is downloaded (called from 'onLoad' in <BODY>)", 'function launch()', '{', 'for (var i = first_image + 1; i <= last_image; i++)', '{', 'theImages[i-first_image] = new Image();', 'theImages[i-first_image].src = modImages[i-first_image];', 'imageNum[i-first_image] = true;', 'document.animation.src = theImages[i-first_image].src;', 'document.control_form.frame_nr.value = i;', '}', '', '// this needs to be done to set the right mode when the page is manually reloaded', 'change_mode (1);', 'fwd();', '}', '', '//===> Check selection status of image in animation loop', 'function checkImage(status,i)', '{', 'if (status == true)', 'imageNum[i] = false;', 'else imageNum[i] = true;', '}', '', '//==> Empty function - used to deal with image buttons rather than HTML buttons', 'function func()', '{', '}', '', '//===> Sets up interface - this is the one function called from the HTML body', 'function animation()', '{', 'count = first_image;', '}', '', '// -->', '</SCRIPT>', '</HEAD>', '', '<BODY BGCOLOR="#003366" onLoad="launch()">', '', '<P ALIGN=LEFT>', '<BR CLEAR=ALL>', '<UL>', '', '</P>', '', '<!-- <CENTER> <P ALIGN=CENTER> -->', '<TABLE ALIGN=left BORDER=2 CELLPADDING=3 CELLSPACING=3>', '<TR>', '<TH ALIGN=left BGCOLOR="#B0C4DE" WIDTH=150> <p class="control">Frame Controls:</p></TH>', '<TH ALIGN=left BGCOLOR="#B0C4DE">&#160;</TH>', '</TR>', '<TR VALIGN="left">', '<TD BGCOLOR="#B0C4DE">', "<IMG SRC='bar_cld.png' width=150 height=600><br>", '<A HREF="JavaScript: func()" onClick="decrementImage(--current_image)"><IMG BORDER=0 ALT="-1"></A>', '<A HREF="JavaScript: func()" onClick="change_mode(1);rev()"><IMG BORDER=0 ALT="Rev"></A>', '<A HREF="JavaScript: func()" onClick="stop()"><IMG BORDER=0 ALT="Stop"></A>', '<A HREF="JavaScript: func()" onClick="change_mode(1);fwd()"><IMG BORDER=0 ALT="Fwd"></A>', '<A HREF="JavaScript: func()" onClick="incrementImage(++current_image)"><IMG BORDER=0 ALT="+1"></A>', '<BR>&#160;<BR>', '<p class="control2">Adjust Speed:<BR>', '<A HREF="JavaScript: func()" onClick="change_speed(delay_step)"><IMG BORDER=0 ALT="slow"></A>', '<A HREF="JavaScript: func()" onClick="change_speed(-delay_step)"><IMG BORDER=0 ALT="fast"></A>', '<BR>&#160;</p>', '<FORM METHOD="POST" NAME="control_form">', '<p class="control2">Frame:', '<INPUT TYPE="text" NAME="frame_nr" VALUE=9 SIZE="2" onFocus="this.select()" onChange="go2image(this.value)"></INPUT><BR>&#160;</p>', '</FORM>', '</TD>', '<TD BGCOLOR="#B0C4DE" ALIGN=left VALIGN=left><IMG NAME="animation" BORDER=0 WIDTH=800 HEIGHT=800 ALT="image"></TD>', '</TR>', '</TABLE>', '</P> </CENTER>', '<BR CLEAR="ALL">', '</BODY>', '</HTML>']
    return template_lines
