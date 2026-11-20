// LaTeX article theme for Markdeep 2.
//
// A faithful port of the Markdeep 1.20 "latex.css" template: a narrow,
// centered, justified Palatino column with LaTeX-like headings, numbered
// sections, italic captions, and the Markdeep 1.20 contents box.
//
// Put this file beside markdeep2.min.js and make it the first line of a
// .md.html document; nothing else is needed. It must come before <plaintext>,
// which is where Markdeep 2's document source begins:
//
//   <meta charset="utf-8"><script src="latex-theme-v2.js" markdeep="2"></script><plaintext>

// The whole port lives in one stylesheet so that the theme stays a single
// file. A data: URL is a normal stylesheet URL to Markdeep's stylesheets map.
const LATEX_CSS = String.raw`
/* -------------------------------------------------------------------------
   Page: Markdeep 1.20's body rule followed by latex.css's overrides.
   ------------------------------------------------------------------------- */
body.md {
    display: block;
    max-width: 75%;
    margin: 0 auto;
    padding: 20px;
    color: #000;
    background-color: #fff;
    font-family: Palatino, Georgia, "Times New Roman", serif;
    font-size: 14px;
    font-weight: 400;
    line-height: 140%;
    text-align: justify;
    -webkit-font-smoothing: antialiased;
    -moz-osx-font-smoothing: grayscale;
    font-smoothing: antialiased;

    /* Markdeep 2 measures the column with these; keep them in step with the
       body so captions and code blocks see the same width as Markdeep 1. */
    --article-padding: 0px;
    --content-max-width: 100%;
}

/* Markdeep 2 lays documents out in a grid of page panes. The 1.20 template is
   one plain column: flatten the panes and remove their padding. */
#main, #content { display: block; }
#article {
    width: auto;
    max-width: none;
    margin: 0;
    padding: 0;
    flex: none;
}
#main_header, #main_footer { padding: 0; }

/* Markdeep 2 puts every code fence and wide table in a scroll viewport. */
.md .scrollable {
    width: auto;
    margin-inline: 0;
    overflow: visible;
}
.md .scrollable > .scrollbar { display: none; }

/* -------------------------------------------------------------------------
   Titles and headings
   latex.css: title 24px normal, padding 170/0/12; h1 18px, h2 15px,
   h3-h6 14px, all Palatino with no rules.
   ------------------------------------------------------------------------- */
.md h1, .md div.title {
    font-family: Palatino, Georgia, "Times New Roman", serif;
    font-size: 24px;
    font-weight: normal;
    line-height: 120%;
    text-align: center;
    color: #000;
    margin: 13.4px 0;
    padding: 170px 0 12px;
    border: none;
    clear: both;
}

.md div.subtitle {
    text-align: center;
    font-size: 15px;
    padding-bottom: 12px;
}

.md h2, .md h3, .md h4, .md h5, .md h6,
.md .nonumberh1, .md .nonumberh2, .md .nonumberh3,
.md .nonumberh4, .md .nonumberh5, .md .nonumberh6 {
    font-family: Palatino, Georgia, "Times New Roman", serif;
    font-weight: bold;
    line-height: inherit;
    text-align: inherit;
    color: #000;
    margin: 13.4px 0;
    padding: 15px 0 3px;
    border: none;
    clear: both;
}

/* Markdeep 1.20 started every section with an empty inline-block anchor,
   which kept the heading's top margin inside the section. Markdeep 2 nests
   the heading in a <section> and the margins collapse, so restore the gap. */
.md section.h1-section, .md section.h2-section, .md section.h3-section,
.md section.h4-section, .md section.h5-section, .md section.h6-section {
    padding-top: 13.4px;
}
.md section.h1-section { padding-top: 0; }
.md section > h1, .md section > h2, .md section > h3,
.md section > h4, .md section > h5, .md section > h6 {
    margin-top: 0;
}
.md section.h1-section > h1 { margin-top: 13.4px; }

.md div.afterTitles { height: 10px; }

.md h2, .md .nonumberh1, .md .nonumberh2 { font-size: 18px; }
.md h3, .md .nonumberh3 { font-size: 15px; }
.md h4, .md h5, .md h6,
.md .nonumberh4, .md .nonumberh5, .md .nonumberh6 { font-size: 14px; }

/* -------------------------------------------------------------------------
   Inline text
   ------------------------------------------------------------------------- */
.md a { font-family: Georgia, Palatino, "Times New Roman", serif; }
.md a:link, .md a:visited { color: #38A; text-decoration: none; }
.md a:link:hover { text-decoration: underline; }

/* latex.css */
.md sup { vertical-align: baseline; position: relative; top: -0.4em; }
.md sub { top: 0.4em; }

/* Markdeep 2 writes *bold* as em.asterisk and __bold__ as strong.underscore;
   the 1.20 stylesheet used a full 700 weight for bold. */
.md em.asterisk { font-style: normal; font-weight: 700; }
.md strong.underscore { font-weight: 700; }
.md strong.asterisk { font-style: italic; font-weight: inherit; }

/* -------------------------------------------------------------------------
   Code
   ------------------------------------------------------------------------- */
.md code, .md pre {
    font-family: Menlo, Consolas, monospace;
    text-align: left;
    line-height: normal;
    white-space: pre;
}

/* latex.css */
.md code { font-size: 82%; }

/* Markdeep 1.20 had no inline-code box. */
.md :not(pre) > code {
    background: none;
    border-radius: 0;
    padding: 0;
}

.md pre.listing {
    width: 100%;
    min-width: 0;
    max-width: 100%;
    box-sizing: content-box;
    position: static;
    left: auto;
    translate: none;
    tab-size: 4;
    -moz-tab-size: 4;
    -o-tab-size: 4;
    counter-reset: line;
    overflow-x: auto;
    resize: horizontal;
}

/* Markdeep 1.20 rules for ~~~ fences. */
.md pre.listing.tilde {
    border: none;
    border-top: 1px solid #CCC;
    border-bottom: 1px solid #CCC;
    border-radius: 0;
    padding: 5px 0 5px 20px;
    margin: 0;
    background: #FCFCFC;
}

/* Markdeep 1.20's highlight.js theme. */
.hljs { display: block; overflow-x: auto; padding: 0.5em; background: #fff; color: #000; -webkit-text-size-adjust: none; }

/* Remove the token rules that only Markdeep 2 adds, so that unstyled token
   classes fall back to the 1.20 defaults. */
.md .hljs-attr, .md .hljs-attribute, .md .hljs-selector-tag, .md .hljs-selector-id,
.md .hljs-selector-class, .md .hljs-selector-attr, .md .hljs-selector-pseudo,
.md .hljs-meta, .md .hljs-section, .md .hljs-bullet, .md .hljs-quote,
.md .hljs-code, .md .hljs-link, .md .hljs-title, .md .hljs-class,
.md .hljs-property { color: inherit; }
.md .hljs-section, .md .hljs-strong { font-weight: inherit; }
.md .hljs-emphasis { font-style: inherit; }
.md .hljs-link { text-decoration: none; }

.hljs-comment { color: #006a00; }
.hljs-keyword { color: #02E; }
.hljs-literal, .nginx .hljs-title { color: #aa0d91; }
.method, .hljs-list .hljs-title, .hljs-tag .hljs-title, .setting .hljs-value,
.hljs-winutils, .tex .hljs-command, .http .hljs-title, .hljs-request,
.hljs-status, .hljs-name { color: #008; }
.hljs-envvar, .tex .hljs-special { color: #660; }
.hljs-string { color: #c41a16; }
.hljs-tag .hljs-value, .hljs-cdata, .hljs-filter .hljs-argument,
.hljs-attr_selector, .apache .hljs-cbracket, .hljs-date, .hljs-regexp { color: #080; }
.hljs-sub .hljs-identifier, .hljs-pi, .hljs-tag, .hljs-tag .hljs-keyword,
.hljs-decorator, .ini .hljs-title, .hljs-shebang, .hljs-prompt,
.hljs-hexcolor, .hljs-rule .hljs-value, .hljs-symbol, .hljs-symbol .hljs-string,
.hljs-number, .css .hljs-function, .hljs-function .hljs-title,
.coffeescript .hljs-attribute { color: #A0C; }
.hljs-function .hljs-title { font-weight: bold; color: #000; }
.hljs-class .hljs-title, .smalltalk .hljs-class, .hljs-type, .hljs-typename,
.hljs-tag .hljs-attribute, .hljs-doctype, .hljs-class .hljs-id, .hljs-built_in,
.setting, .hljs-params, .clojure .hljs-attribute { color: #5c2699; }
.hljs-variable { color: #3f6e74; }
.css .hljs-tag, .hljs-rule .hljs-property, .hljs-pseudo, .hljs-subst { color: #000; }
.css .hljs-class, .css .hljs-id { color: #9b703f; }
.hljs-value .hljs-important { color: #ff7700; font-weight: bold; }
.hljs-rule .hljs-keyword { color: #c5af75; }
.hljs-annotation, .apache .hljs-sqbracket, .nginx .hljs-built_in { color: #9b859d; }
.hljs-preprocessor, .hljs-preprocessor *, .hljs-pragma { color: #643820; }
.tex .hljs-formula { background-color: #eee; font-style: italic; }
.diff .hljs-header, .hljs-chunk { color: #808080; font-weight: bold; }
.diff .hljs-change { background-color: #bccff9; }
.hljs-addition { background-color: #baeeba; }
.hljs-deletion { background-color: #ffc8bd; }
.hljs-comment .hljs-doctag { font-weight: bold; }
.method .hljs-id { color: #000; }

/* -------------------------------------------------------------------------
   Lists, quotes, rules
   ------------------------------------------------------------------------- */
.md ul li.asterisk:before, .md ul li.minus:before {
    content: "\25CF";
    display: block;
    float: left;
    width: 1.2em;
    margin-left: -18px;
}
.md ul li.asterisk, .md ul li.minus {
    list-style-type: none;
    margin-bottom: 10px;
    padding-left: 4pt;
}
.md ol li {
    margin-bottom: 10px;
    padding-left: 4pt;
}

/* Markdeep 2 names task-list items checked_box/unchecked_box. */
.md li.checked_box, .md li.unchecked_box { opacity: 1; overflow: visible; text-indent: -1.2em; }
.md ul > li.checked_box, .md ul > li.unchecked_box { list-style: none; }
.md ol > li.checked_box, .md ol > li.unchecked_box { padding-left: 0; }
.md li.checked_box:before, .md li.unchecked_box:before {
    content: "\2611";
    display: block;
    float: left;
    width: 1em;
    margin-left: 0;
    font-size: 120%;
    font-weight: inherit;
    line-height: inherit;
    position: static;
    top: auto;
}
.md li.unchecked_box:before { content: "\2610"; }

.md hr {
    height: 0;
    border: none;
    border-top: 0.1pt solid #777;
}

.md .admonition.tip::before { color: rgba(26,128,46,.8); }
.md .admonition.warn::before, .md .admonition.warning::before { color: rgba(128,73,0,.8); }
.md .admonition.error::before { color: rgba(128,12,34,.8); }

/* -------------------------------------------------------------------------
   Tables, captions, figures
   ------------------------------------------------------------------------- */
.md table.table th {
    color: #FFF;
    background-color: #AAA;
    border: 1px solid #888;
    padding: 8px 15px;
}
.md table.table td {
    padding: 5px 15px;
    border: 1px solid #888;
}
.md table.table tr:nth-child(even) { background: #EEE; }

/* latex.css targeted only div captions; rendered figures use a span. */
.md div.imagecaption { font-style: normal; }
.md center.wide .imagecaption { max-width: 100%; }

.md svg.diagram { zoom: 80%; stroke-width: 2px; }

/* -------------------------------------------------------------------------
   Notes, bibliography, footer
   ------------------------------------------------------------------------- */
.md .endnote { font-size: 11px; }
.md .bib b { font-weight: normal; }

.markdeepFooter {
    font-size: 9px;
    text-align: right;
    padding-top: 80px;
    color: #999;
    clear: both;
    font-style: italic;
}

/* -------------------------------------------------------------------------
   Contents, rebuilt by markdeepLatexContents() in the Markdeep 1.20 format
   ------------------------------------------------------------------------- */
.md .latexTOCTop { display: none; }
.md .latexTOCMedium {
    float: right;
    font-family: Verdana, Helvetica, Arial, sans-serif;
    font-size: 12px;
    line-height: 15px;
    border-left: 1px solid #CCC;
    padding: 15px 0 3px 15px;
    margin: 15px 0 15px 25px;
    clear: both;
}
.md .latexTOCMedium .level1 { font-weight: 600; }
.md .latexTOCLong .level1 {
    font-weight: 600;
    display: block;
    padding-top: 12px;
    margin: 0 0 -20px;
}
.md .latexTOCShort {
    font-family: Verdana, Helvetica, Arial, sans-serif;
    text-align: center;
    font-weight: bold;
    margin: 15px 13.4px 13.4px;
    padding: 15px 0 3px;
    font-size: 14px;
    clear: both;
}
.md .latexTOCHeader {
    font-family: Verdana, Helvetica, Arial, sans-serif;
    font-size: 20px;
    font-weight: bold;
    line-height: inherit;
    text-align: inherit;
    color: #000;
    margin: 13.4px 0;s
    padding: 15px 0 3px;
    border-top: none;
    border-bottom: 3px solid;
    clear: both;
}

/* -------------------------------------------------------------------------
   Print (latex.css)
   ------------------------------------------------------------------------- */
@media print {
    body.md { max-width: 64%; }
    .md h1, .md div.title { padding-top: 150px; }
}
`;

// Build the Markdeep 1.20 contents box after Markdeep 2 has rendered.
const markdeepLatexContents = function () {
    const article = document.getElementById('article');
    if (!article || article.querySelector('.latexTOC')) return;

    // Markdeep 1.20 kept an empty paragraph between the title and the first
    // section, worth 1em. Markdeep 2 drops it, so put the gap back.
    const addTitleGap = function () {
        const afterTitles = article.querySelector('div.afterTitles');
        let next = afterTitles ? afterTitles.nextElementSibling : null;
        while (next && (next.classList.contains('decoration') ||
            (next.tagName === 'P' && next.textContent.trim() === ''))) {
            next = next.nextElementSibling;
        }
        if (next) next.style.marginTop = '1em';
    };

    const headings = Array.from(article.querySelectorAll('h2, h3, h4, h5, h6'));
    if (headings.length === 0) {
        addTitleGap();
        return;
    }

    // Markdeep 1.20's automatic choice: none, short, medium, or long.
    // Levels are shifted by one in Markdeep 2 because h1 is the title.
    const html = article.innerHTML;
    const total = headings.length;
    const topLevel = headings.filter(function (h) { return h.tagName === 'H2'; }).length;
    const firstHeading = html.search(/<h[2-6][\s>]/);
    const firstAt = firstHeading < 0 ? 0 : firstHeading;
    let tocStyle;
    if ((total < 4 && topLevel <= 1) || html.length < 2048) tocStyle = 'none';
    else if (topLevel < 7 && total / topLevel < 2.5) tocStyle = 'short';
    else if (firstHeading < 0 || firstAt / 55 > total) tocStyle = 'medium';
    else tocStyle = 'long';
    if (tocStyle === 'none') {
        addTitleGap();
        return;
    }

    let entries = '';
    let shortTOC = '';
    headings.forEach(function (heading) {
        const level = parseInt(heading.tagName.charAt(1), 10) - 1;
        if (level > 3) return;                  // Markdeep 1.20's tocDepth 3

        const target = heading.querySelector('a.target');
        const href = '#' + (target ? target.getAttribute('name') : '');
        const copy = heading.cloneNode(true);
        copy.querySelectorAll('a.target').forEach(function (a) { a.remove(); });
        const title = copy.innerHTML.trim();
        const number = heading.getAttribute('data-section-number');

        entries += '&nbsp;&nbsp;'.repeat(level - 1) +
            '<a href="' + href + '" target="_self" class="level' + level + '">' +
            (number ? '<span class="tocNumber">' + number + '&nbsp; </span>' : '') +
            title + '</a><br>\n';
        if (level === 1) {
            shortTOC += (shortTOC ? ' &middot; ' : '') +
                '<a href="' + href + '" target="_self">' + title + '</a>';
        }
    });

    const toc = document.createElement('div');
    if (tocStyle === 'short') {
        toc.className = 'latexTOC latexTOCShort';
        toc.innerHTML = shortTOC;
    } else {
        const medium = tocStyle === 'medium';
        toc.className = 'latexTOC ' + (medium ? 'latexTOCMedium' : 'latexTOCLong');
        toc.innerHTML = (medium ? '<center><b>Table of Contents</b></center>' : '<div class="latexTOCHeader">Table of Contents</div>') +
            '<p><a href="#" class="latexTOCTop" target="_self">(Top)</a><br>\n' + entries + '</p>';
    }

    if (tocStyle === 'long') {
        // Markdeep 1.20 put the long contents block in the flow just before the
        // first section. Insert before the section wrapper so its top padding
        // does not push the contents down.
        const section = headings[0].closest('section') || headings[0];
        section.parentNode.insertBefore(toc, section);
    } else {
        const afterTitles = article.querySelector('div.afterTitles');
        const title = article.querySelector('h1');
        const anchor = afterTitles || title;
        if (anchor) anchor.parentNode.insertBefore(toc, anchor.nextSibling);
        else article.insertBefore(toc, article.firstChild);
    }
};

( window.markdeepOptions ||=[] ).push({
    stylesheets: {latex: 'data:text/css;charset=utf-8,' + encodeURIComponent(LATEX_CSS)},

    // Markdeep 1.20 defaults that affect layout rather than type.
    captionAbove: {diagram: false, image: false, table: false, listing: false},
    tocDepth: 4,                 // 1.20's tocDepth 3, shifted one heading level
    showSectionNumbers: true,
    stripSectionNumbers: false,
    showLinkPreviews: false,

    // Markdeep 1.20 had one light color scheme.
    initialColorScheme: 'light',
    enableColorSchemeChoice: false,

    // The 1.20 contents box is rebuilt in the text column above.
    main_nav: {align: 'none'},

    postRenderCallback: markdeepLatexContents,
});

// Mandatory activation line. DO NOT MODIFY OR RELOCATE
((d,o)=>d.head.appendChild(((m)=>{m.dataset.theme=d.currentScript.src;return m;})(Object.assign(d.createElement('script'),{src:o.find((e)=>e.markdeepPath)?.markdeepPath??d.currentScript.src.replace(/[^/]*$/,'markdeep2.min.js')}))))(document,window.markdeepOptions);
