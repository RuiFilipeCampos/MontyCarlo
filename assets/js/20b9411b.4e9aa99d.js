"use strict";(self.webpackChunkdev=self.webpackChunkdev||[]).push([[722],{3905:function(a,e,t){t.d(e,{Zo:function(){return i},kt:function(){return k}});var s=t(7294);function n(a,e,t){return e in a?Object.defineProperty(a,e,{value:t,enumerable:!0,configurable:!0,writable:!0}):a[e]=t,a}function m(a,e){var t=Object.keys(a);if(Object.getOwnPropertySymbols){var s=Object.getOwnPropertySymbols(a);e&&(s=s.filter((function(e){return Object.getOwnPropertyDescriptor(a,e).enumerable}))),t.push.apply(t,s)}return t}function p(a){for(var e=1;e<arguments.length;e++){var t=null!=arguments[e]?arguments[e]:{};e%2?m(Object(t),!0).forEach((function(e){n(a,e,t[e])})):Object.getOwnPropertyDescriptors?Object.defineProperties(a,Object.getOwnPropertyDescriptors(t)):m(Object(t)).forEach((function(e){Object.defineProperty(a,e,Object.getOwnPropertyDescriptor(t,e))}))}return a}function r(a,e){if(null==a)return{};var t,s,n=function(a,e){if(null==a)return{};var t,s,n={},m=Object.keys(a);for(s=0;s<m.length;s++)t=m[s],e.indexOf(t)>=0||(n[t]=a[t]);return n}(a,e);if(Object.getOwnPropertySymbols){var m=Object.getOwnPropertySymbols(a);for(s=0;s<m.length;s++)t=m[s],e.indexOf(t)>=0||Object.prototype.propertyIsEnumerable.call(a,t)&&(n[t]=a[t])}return n}var l=s.createContext({}),N=function(a){var e=s.useContext(l),t=e;return a&&(t="function"==typeof a?a(e):p(p({},e),a)),t},i=function(a){var e=N(a.components);return s.createElement(l.Provider,{value:e},a.children)},c={inlineCode:"code",wrapper:function(a){var e=a.children;return s.createElement(s.Fragment,{},e)}},o=s.forwardRef((function(a,e){var t=a.components,n=a.mdxType,m=a.originalType,l=a.parentName,i=r(a,["components","mdxType","originalType","parentName"]),o=N(t),k=n,h=o["".concat(l,".").concat(k)]||o[k]||c[k]||m;return t?s.createElement(h,p(p({ref:e},i),{},{components:t})):s.createElement(h,p({ref:e},i))}));function k(a,e){var t=arguments,n=e&&e.mdxType;if("string"==typeof a||n){var m=t.length,p=new Array(m);p[0]=o;var r={};for(var l in e)hasOwnProperty.call(e,l)&&(r[l]=e[l]);r.originalType=a,r.mdxType="string"==typeof a?a:n,p[1]=r;for(var N=2;N<m;N++)p[N]=t[N];return s.createElement.apply(null,p)}return s.createElement.apply(null,t)}o.displayName="MDXCreateElement"},4744:function(a,e,t){t.r(e),t.d(e,{frontMatter:function(){return r},contentTitle:function(){return l},metadata:function(){return N},toc:function(){return i},default:function(){return o}});var s=t(7462),n=t(3366),m=(t(7294),t(3905)),p=["components"],r={},l="Coherent Scattering",N={unversionedId:"photons/coherent",id:"photons/coherent",isDocsHomePage:!1,title:"Coherent Scattering",description:"Thomson Scattering",source:"@site/physics/photons/coherent.md",sourceDirName:"photons",slug:"/photons/coherent",permalink:"/MontyCarlo/physics/photons/coherent",tags:[],version:"current",frontMatter:{},sidebar:"tutorialSidebar",previous:{title:"Methods for Random Sampling",permalink:"/MontyCarlo/physics/random_sampling/uniform"},next:{title:"Compton Scattering",permalink:"/MontyCarlo/physics/photons/incoherent"}},i=[{value:"Thomson Scattering",id:"thomson-scattering",children:[],level:2},{value:"Rayleigh Scattering",id:"rayleigh-scattering",children:[],level:2}],c={toc:i};function o(a){var e=a.components,t=(0,n.Z)(a,p);return(0,m.kt)("wrapper",(0,s.Z)({},c,t,{components:e,mdxType:"MDXLayout"}),(0,m.kt)("h1",{id:"coherent-scattering"},"Coherent Scattering"),(0,m.kt)("h2",{id:"thomson-scattering"},"Thomson Scattering"),(0,m.kt)("p",null,"Consider an incident electromagnetic plane wave on a free charged particle that is initially at rest. The EM wave can be written in versor notation"),(0,m.kt)("div",{className:"math math-display"},(0,m.kt)("span",{parentName:"div",className:"katex-display"},(0,m.kt)("span",{parentName:"span",className:"katex"},(0,m.kt)("span",{parentName:"span",className:"katex-mathml"},(0,m.kt)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML",display:"block"},(0,m.kt)("semantics",{parentName:"math"},(0,m.kt)("mrow",{parentName:"semantics"},(0,m.kt)("mover",{parentName:"mrow",accent:"true"},(0,m.kt)("mi",{parentName:"mover",mathvariant:"bold"},"E"),(0,m.kt)("mo",{parentName:"mover"},"~")),(0,m.kt)("mo",{parentName:"mrow"},"="),(0,m.kt)("mover",{parentName:"mrow",accent:"true"},(0,m.kt)("mi",{parentName:"mover",mathvariant:"bold"},"z"),(0,m.kt)("mo",{parentName:"mover"},"^")),(0,m.kt)("msub",{parentName:"mrow"},(0,m.kt)("mi",{parentName:"msub"},"E"),(0,m.kt)("mn",{parentName:"msub"},"0")),(0,m.kt)("mi",{parentName:"mrow"},"exp"),(0,m.kt)("mo",{parentName:"mrow"},"\u2061"),(0,m.kt)("mrow",{parentName:"mrow"},(0,m.kt)("mo",{parentName:"mrow",fence:"true"},"("),(0,m.kt)("mi",{parentName:"mrow"},"i"),(0,m.kt)("mo",{parentName:"mrow",stretchy:"false"},"("),(0,m.kt)("mi",{parentName:"mrow",mathvariant:"bold"},"k"),(0,m.kt)("mo",{parentName:"mrow"},"\u22c5"),(0,m.kt)("mi",{parentName:"mrow",mathvariant:"bold"},"r"),(0,m.kt)("mo",{parentName:"mrow"},"\u2212"),(0,m.kt)("mi",{parentName:"mrow"},"w"),(0,m.kt)("mi",{parentName:"mrow"},"t"),(0,m.kt)("mo",{parentName:"mrow",stretchy:"false"},")"),(0,m.kt)("mo",{parentName:"mrow",fence:"true"},")"))),(0,m.kt)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"\\tilde \\mathbf{E} = \\hat \\mathbf{ z} E_0 \\exp \\left ( i (\\mathbf{k} \\cdot \\mathbf{r} - wt) \\right )")))),(0,m.kt)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},(0,m.kt)("span",{parentName:"span",className:"base"},(0,m.kt)("span",{parentName:"span",className:"strut",style:{height:"0.9229699999999998em",verticalAlign:"0em"}}),(0,m.kt)("span",{parentName:"span",className:"mord accent"},(0,m.kt)("span",{parentName:"span",className:"vlist-t"},(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.9229699999999998em"}},(0,m.kt)("span",{parentName:"span",style:{top:"-3em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathbf"},"E"))),(0,m.kt)("span",{parentName:"span",style:{top:"-3.6051100000000003em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"accent-body",style:{left:"-0.16666em"}},(0,m.kt)("span",{parentName:"span",className:"mord"},"~"))))))),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),(0,m.kt)("span",{parentName:"span",className:"mrel"},"="),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),(0,m.kt)("span",{parentName:"span",className:"base"},(0,m.kt)("span",{parentName:"span",className:"strut",style:{height:"1em",verticalAlign:"-0.25em"}}),(0,m.kt)("span",{parentName:"span",className:"mord accent"},(0,m.kt)("span",{parentName:"span",className:"vlist-t"},(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.70788em"}},(0,m.kt)("span",{parentName:"span",style:{top:"-3em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathbf"},"z"))),(0,m.kt)("span",{parentName:"span",style:{top:"-3.01344em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"accent-body",style:{left:"-0.19444em"}},(0,m.kt)("span",{parentName:"span",className:"mord"},"^"))))))),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathnormal",style:{marginRight:"0.05764em"}},"E"),(0,m.kt)("span",{parentName:"span",className:"msupsub"},(0,m.kt)("span",{parentName:"span",className:"vlist-t vlist-t2"},(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.30110799999999993em"}},(0,m.kt)("span",{parentName:"span",style:{top:"-2.5500000000000003em",marginLeft:"-0.05764em",marginRight:"0.05em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),(0,m.kt)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},(0,m.kt)("span",{parentName:"span",className:"mord mtight"},"0")))),(0,m.kt)("span",{parentName:"span",className:"vlist-s"},"\u200b")),(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.15em"}},(0,m.kt)("span",{parentName:"span"})))))),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),(0,m.kt)("span",{parentName:"span",className:"mop"},"exp"),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),(0,m.kt)("span",{parentName:"span",className:"minner"},(0,m.kt)("span",{parentName:"span",className:"mopen delimcenter",style:{top:"0em"}},"("),(0,m.kt)("span",{parentName:"span",className:"mord mathnormal"},"i"),(0,m.kt)("span",{parentName:"span",className:"mopen"},"("),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathbf"},"k")),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}}),(0,m.kt)("span",{parentName:"span",className:"mbin"},"\u22c5"),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathbf"},"r")),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}}),(0,m.kt)("span",{parentName:"span",className:"mbin"},"\u2212"),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}}),(0,m.kt)("span",{parentName:"span",className:"mord mathnormal",style:{marginRight:"0.02691em"}},"w"),(0,m.kt)("span",{parentName:"span",className:"mord mathnormal"},"t"),(0,m.kt)("span",{parentName:"span",className:"mclose"},")"),(0,m.kt)("span",{parentName:"span",className:"mclose delimcenter",style:{top:"0em"}},")"))))))),(0,m.kt)("p",null,"where the actual real value electric field can be extracted by taking the real part of ",(0,m.kt)("span",{parentName:"p",className:"math math-inline"},(0,m.kt)("span",{parentName:"span",className:"katex"},(0,m.kt)("span",{parentName:"span",className:"katex-mathml"},(0,m.kt)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML"},(0,m.kt)("semantics",{parentName:"math"},(0,m.kt)("mrow",{parentName:"semantics"},(0,m.kt)("mover",{parentName:"mrow",accent:"true"},(0,m.kt)("mi",{parentName:"mover",mathvariant:"bold"},"E"),(0,m.kt)("mo",{parentName:"mover"},"~"))),(0,m.kt)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"\\tilde \\mathbf{E}")))),(0,m.kt)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},(0,m.kt)("span",{parentName:"span",className:"base"},(0,m.kt)("span",{parentName:"span",className:"strut",style:{height:"0.9229699999999998em",verticalAlign:"0em"}}),(0,m.kt)("span",{parentName:"span",className:"mord accent"},(0,m.kt)("span",{parentName:"span",className:"vlist-t"},(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.9229699999999998em"}},(0,m.kt)("span",{parentName:"span",style:{top:"-3em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathbf"},"E"))),(0,m.kt)("span",{parentName:"span",style:{top:"-3.6051100000000003em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"accent-body",style:{left:"-0.16666em"}},(0,m.kt)("span",{parentName:"span",className:"mord"},"~"))))))))))),", that is,"),(0,m.kt)("div",{className:"math math-display"},(0,m.kt)("span",{parentName:"div",className:"katex-display"},(0,m.kt)("span",{parentName:"span",className:"katex"},(0,m.kt)("span",{parentName:"span",className:"katex-mathml"},(0,m.kt)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML",display:"block"},(0,m.kt)("semantics",{parentName:"math"},(0,m.kt)("mrow",{parentName:"semantics"},(0,m.kt)("mi",{parentName:"mrow",mathvariant:"bold"},"E"),(0,m.kt)("mo",{parentName:"mrow"},"="),(0,m.kt)("mtext",{parentName:"mrow"},"Re"),(0,m.kt)("mo",{parentName:"mrow",stretchy:"false"},"("),(0,m.kt)("mover",{parentName:"mrow",accent:"true"},(0,m.kt)("mi",{parentName:"mover",mathvariant:"bold"},"E"),(0,m.kt)("mo",{parentName:"mover"},"~")),(0,m.kt)("mo",{parentName:"mrow",stretchy:"false"},")"),(0,m.kt)("mi",{parentName:"mrow",mathvariant:"normal"},".")),(0,m.kt)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"\\mathbf{E} = \\textrm{Re}( \\tilde \\mathbf{E} ).")))),(0,m.kt)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},(0,m.kt)("span",{parentName:"span",className:"base"},(0,m.kt)("span",{parentName:"span",className:"strut",style:{height:"0.68611em",verticalAlign:"0em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathbf"},"E")),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),(0,m.kt)("span",{parentName:"span",className:"mrel"},"="),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),(0,m.kt)("span",{parentName:"span",className:"base"},(0,m.kt)("span",{parentName:"span",className:"strut",style:{height:"1.1729699999999998em",verticalAlign:"-0.25em"}}),(0,m.kt)("span",{parentName:"span",className:"mord text"},(0,m.kt)("span",{parentName:"span",className:"mord textrm"},"Re")),(0,m.kt)("span",{parentName:"span",className:"mopen"},"("),(0,m.kt)("span",{parentName:"span",className:"mord accent"},(0,m.kt)("span",{parentName:"span",className:"vlist-t"},(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.9229699999999998em"}},(0,m.kt)("span",{parentName:"span",style:{top:"-3em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathbf"},"E"))),(0,m.kt)("span",{parentName:"span",style:{top:"-3.6051100000000003em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"accent-body",style:{left:"-0.16666em"}},(0,m.kt)("span",{parentName:"span",className:"mord"},"~"))))))),(0,m.kt)("span",{parentName:"span",className:"mclose"},")"),(0,m.kt)("span",{parentName:"span",className:"mord"},".")))))),(0,m.kt)("p",null,"Note that since the field is always pointing in the z-direction, the charge will be set in motion in that direction. Meaning that this is a 1d problem and only the z-coordinate needs to be considered. Using Newtons Law,"),(0,m.kt)("div",{className:"math math-display"},(0,m.kt)("span",{parentName:"div",className:"katex-display"},(0,m.kt)("span",{parentName:"span",className:"katex"},(0,m.kt)("span",{parentName:"span",className:"katex-mathml"},(0,m.kt)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML",display:"block"},(0,m.kt)("semantics",{parentName:"math"},(0,m.kt)("mrow",{parentName:"semantics"},(0,m.kt)("mi",{parentName:"mrow"},"m"),(0,m.kt)("mi",{parentName:"mrow",mathvariant:"bold"},"a"),(0,m.kt)("mo",{parentName:"mrow"},"\u22c5"),(0,m.kt)("mover",{parentName:"mrow",accent:"true"},(0,m.kt)("mi",{parentName:"mover",mathvariant:"bold"},"z"),(0,m.kt)("mo",{parentName:"mover"},"^")),(0,m.kt)("mo",{parentName:"mrow"},"="),(0,m.kt)("mi",{parentName:"mrow"},"e"),(0,m.kt)("mtext",{parentName:"mrow"},"Re"),(0,m.kt)("mo",{parentName:"mrow",stretchy:"false"},"("),(0,m.kt)("mover",{parentName:"mrow",accent:"true"},(0,m.kt)("mi",{parentName:"mover",mathvariant:"bold"},"E"),(0,m.kt)("mo",{parentName:"mover"},"~")),(0,m.kt)("mo",{parentName:"mrow",stretchy:"false"},")"),(0,m.kt)("mi",{parentName:"mrow",mathvariant:"normal"},".")),(0,m.kt)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"m \\mathbf{a} \\cdot \\hat \\mathbf{z} = e \\textrm{Re}( \\tilde \\mathbf{E} ).")))),(0,m.kt)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},(0,m.kt)("span",{parentName:"span",className:"base"},(0,m.kt)("span",{parentName:"span",className:"strut",style:{height:"0.44445em",verticalAlign:"0em"}}),(0,m.kt)("span",{parentName:"span",className:"mord mathnormal"},"m"),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathbf"},"a")),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}}),(0,m.kt)("span",{parentName:"span",className:"mbin"},"\u22c5"),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}})),(0,m.kt)("span",{parentName:"span",className:"base"},(0,m.kt)("span",{parentName:"span",className:"strut",style:{height:"0.70788em",verticalAlign:"0em"}}),(0,m.kt)("span",{parentName:"span",className:"mord accent"},(0,m.kt)("span",{parentName:"span",className:"vlist-t"},(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.70788em"}},(0,m.kt)("span",{parentName:"span",style:{top:"-3em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathbf"},"z"))),(0,m.kt)("span",{parentName:"span",style:{top:"-3.01344em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"accent-body",style:{left:"-0.19444em"}},(0,m.kt)("span",{parentName:"span",className:"mord"},"^"))))))),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),(0,m.kt)("span",{parentName:"span",className:"mrel"},"="),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),(0,m.kt)("span",{parentName:"span",className:"base"},(0,m.kt)("span",{parentName:"span",className:"strut",style:{height:"1.1729699999999998em",verticalAlign:"-0.25em"}}),(0,m.kt)("span",{parentName:"span",className:"mord mathnormal"},"e"),(0,m.kt)("span",{parentName:"span",className:"mord text"},(0,m.kt)("span",{parentName:"span",className:"mord textrm"},"Re")),(0,m.kt)("span",{parentName:"span",className:"mopen"},"("),(0,m.kt)("span",{parentName:"span",className:"mord accent"},(0,m.kt)("span",{parentName:"span",className:"vlist-t"},(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.9229699999999998em"}},(0,m.kt)("span",{parentName:"span",style:{top:"-3em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathbf"},"E"))),(0,m.kt)("span",{parentName:"span",style:{top:"-3.6051100000000003em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"accent-body",style:{left:"-0.16666em"}},(0,m.kt)("span",{parentName:"span",className:"mord"},"~"))))))),(0,m.kt)("span",{parentName:"span",className:"mclose"},")"),(0,m.kt)("span",{parentName:"span",className:"mord"},".")))))),(0,m.kt)("h2",{id:"rayleigh-scattering"},"Rayleigh Scattering"),(0,m.kt)("p",null,"The DCS:"),(0,m.kt)("div",{className:"math math-display"},(0,m.kt)("span",{parentName:"div",className:"katex-display"},(0,m.kt)("span",{parentName:"span",className:"katex"},(0,m.kt)("span",{parentName:"span",className:"katex-mathml"},(0,m.kt)("math",{parentName:"span",xmlns:"http://www.w3.org/1998/Math/MathML",display:"block"},(0,m.kt)("semantics",{parentName:"math"},(0,m.kt)("mrow",{parentName:"semantics"},(0,m.kt)("mfrac",{parentName:"mrow"},(0,m.kt)("mrow",{parentName:"mfrac"},(0,m.kt)("mi",{parentName:"mrow"},"d"),(0,m.kt)("mi",{parentName:"mrow"},"\u03c3")),(0,m.kt)("mrow",{parentName:"mfrac"},(0,m.kt)("mi",{parentName:"mrow"},"d"),(0,m.kt)("mi",{parentName:"mrow",mathvariant:"normal"},"\u03a9"))),(0,m.kt)("mo",{parentName:"mrow"},"="),(0,m.kt)("msubsup",{parentName:"mrow"},(0,m.kt)("mi",{parentName:"msubsup"},"r"),(0,m.kt)("mi",{parentName:"msubsup"},"e"),(0,m.kt)("mn",{parentName:"msubsup"},"2")),(0,m.kt)("mfrac",{parentName:"mrow"},(0,m.kt)("mrow",{parentName:"mfrac"},(0,m.kt)("mn",{parentName:"mrow"},"1"),(0,m.kt)("mo",{parentName:"mrow"},"+"),(0,m.kt)("msup",{parentName:"mrow"},(0,m.kt)("mo",{parentName:"msup"},(0,m.kt)("mi",{parentName:"mo"},"cos"),(0,m.kt)("mo",{parentName:"mo"},"\u2061")),(0,m.kt)("mn",{parentName:"msup"},"2")),(0,m.kt)("mi",{parentName:"mrow"},"\u03b8")),(0,m.kt)("mn",{parentName:"mfrac"},"2")),(0,m.kt)("msub",{parentName:"mrow"},(0,m.kt)("mi",{parentName:"msub"},"F"),(0,m.kt)("mi",{parentName:"msub"},"Z")),(0,m.kt)("mo",{parentName:"mrow",stretchy:"false"},"("),(0,m.kt)("mi",{parentName:"mrow"},"q"),(0,m.kt)("msup",{parentName:"mrow"},(0,m.kt)("mo",{parentName:"msup",stretchy:"false"},")"),(0,m.kt)("mn",{parentName:"msup"},"2")),(0,m.kt)("mo",{parentName:"mrow",separator:"true"},",")),(0,m.kt)("annotation",{parentName:"semantics",encoding:"application/x-tex"},"\\frac{d\\sigma}{d\\Omega} = r_e^2 \\frac{1 + \\cos^2 \\theta}{2} F_Z(q)^2,")))),(0,m.kt)("span",{parentName:"span",className:"katex-html","aria-hidden":"true"},(0,m.kt)("span",{parentName:"span",className:"base"},(0,m.kt)("span",{parentName:"span",className:"strut",style:{height:"2.05744em",verticalAlign:"-0.686em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mopen nulldelimiter"}),(0,m.kt)("span",{parentName:"span",className:"mfrac"},(0,m.kt)("span",{parentName:"span",className:"vlist-t vlist-t2"},(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"1.37144em"}},(0,m.kt)("span",{parentName:"span",style:{top:"-2.314em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathnormal"},"d"),(0,m.kt)("span",{parentName:"span",className:"mord"},"\u03a9"))),(0,m.kt)("span",{parentName:"span",style:{top:"-3.23em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"frac-line",style:{borderBottomWidth:"0.04em"}})),(0,m.kt)("span",{parentName:"span",style:{top:"-3.677em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathnormal"},"d"),(0,m.kt)("span",{parentName:"span",className:"mord mathnormal",style:{marginRight:"0.03588em"}},"\u03c3")))),(0,m.kt)("span",{parentName:"span",className:"vlist-s"},"\u200b")),(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.686em"}},(0,m.kt)("span",{parentName:"span"}))))),(0,m.kt)("span",{parentName:"span",className:"mclose nulldelimiter"})),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}}),(0,m.kt)("span",{parentName:"span",className:"mrel"},"="),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2777777777777778em"}})),(0,m.kt)("span",{parentName:"span",className:"base"},(0,m.kt)("span",{parentName:"span",className:"strut",style:{height:"2.177108em",verticalAlign:"-0.686em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathnormal",style:{marginRight:"0.02778em"}},"r"),(0,m.kt)("span",{parentName:"span",className:"msupsub"},(0,m.kt)("span",{parentName:"span",className:"vlist-t vlist-t2"},(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.8641079999999999em"}},(0,m.kt)("span",{parentName:"span",style:{top:"-2.4530000000000003em",marginLeft:"-0.02778em",marginRight:"0.05em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),(0,m.kt)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},(0,m.kt)("span",{parentName:"span",className:"mord mathnormal mtight"},"e"))),(0,m.kt)("span",{parentName:"span",style:{top:"-3.113em",marginRight:"0.05em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),(0,m.kt)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},(0,m.kt)("span",{parentName:"span",className:"mord mtight"},"2")))),(0,m.kt)("span",{parentName:"span",className:"vlist-s"},"\u200b")),(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.247em"}},(0,m.kt)("span",{parentName:"span"})))))),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mopen nulldelimiter"}),(0,m.kt)("span",{parentName:"span",className:"mfrac"},(0,m.kt)("span",{parentName:"span",className:"vlist-t vlist-t2"},(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"1.491108em"}},(0,m.kt)("span",{parentName:"span",style:{top:"-2.314em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord"},"2"))),(0,m.kt)("span",{parentName:"span",style:{top:"-3.23em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"frac-line",style:{borderBottomWidth:"0.04em"}})),(0,m.kt)("span",{parentName:"span",style:{top:"-3.677em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"3em"}}),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord"},"1"),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}}),(0,m.kt)("span",{parentName:"span",className:"mbin"},"+"),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.2222222222222222em"}}),(0,m.kt)("span",{parentName:"span",className:"mop"},(0,m.kt)("span",{parentName:"span",className:"mop"},"cos"),(0,m.kt)("span",{parentName:"span",className:"msupsub"},(0,m.kt)("span",{parentName:"span",className:"vlist-t"},(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.8141079999999999em"}},(0,m.kt)("span",{parentName:"span",style:{top:"-3.063em",marginRight:"0.05em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),(0,m.kt)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},(0,m.kt)("span",{parentName:"span",className:"mord mtight"},"2")))))))),(0,m.kt)("span",{parentName:"span",className:"mspace",style:{marginRight:"0.16666666666666666em"}}),(0,m.kt)("span",{parentName:"span",className:"mord mathnormal",style:{marginRight:"0.02778em"}},"\u03b8")))),(0,m.kt)("span",{parentName:"span",className:"vlist-s"},"\u200b")),(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.686em"}},(0,m.kt)("span",{parentName:"span"}))))),(0,m.kt)("span",{parentName:"span",className:"mclose nulldelimiter"})),(0,m.kt)("span",{parentName:"span",className:"mord"},(0,m.kt)("span",{parentName:"span",className:"mord mathnormal",style:{marginRight:"0.13889em"}},"F"),(0,m.kt)("span",{parentName:"span",className:"msupsub"},(0,m.kt)("span",{parentName:"span",className:"vlist-t vlist-t2"},(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.32833099999999993em"}},(0,m.kt)("span",{parentName:"span",style:{top:"-2.5500000000000003em",marginLeft:"-0.13889em",marginRight:"0.05em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),(0,m.kt)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},(0,m.kt)("span",{parentName:"span",className:"mord mathnormal mtight",style:{marginRight:"0.07153em"}},"Z")))),(0,m.kt)("span",{parentName:"span",className:"vlist-s"},"\u200b")),(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.15em"}},(0,m.kt)("span",{parentName:"span"})))))),(0,m.kt)("span",{parentName:"span",className:"mopen"},"("),(0,m.kt)("span",{parentName:"span",className:"mord mathnormal",style:{marginRight:"0.03588em"}},"q"),(0,m.kt)("span",{parentName:"span",className:"mclose"},(0,m.kt)("span",{parentName:"span",className:"mclose"},")"),(0,m.kt)("span",{parentName:"span",className:"msupsub"},(0,m.kt)("span",{parentName:"span",className:"vlist-t"},(0,m.kt)("span",{parentName:"span",className:"vlist-r"},(0,m.kt)("span",{parentName:"span",className:"vlist",style:{height:"0.8641079999999999em"}},(0,m.kt)("span",{parentName:"span",style:{top:"-3.113em",marginRight:"0.05em"}},(0,m.kt)("span",{parentName:"span",className:"pstrut",style:{height:"2.7em"}}),(0,m.kt)("span",{parentName:"span",className:"sizing reset-size6 size3 mtight"},(0,m.kt)("span",{parentName:"span",className:"mord mtight"},"2")))))))),(0,m.kt)("span",{parentName:"span",className:"mpunct"},",")))))))}o.isMDXComponent=!0}}]);