// @ts-check
// Note: type annotations allow type checking and IDEs autocompletion

const lightCodeTheme = require('prism-react-renderer/themes/github');
const darkCodeTheme = require('prism-react-renderer/themes/dracula');
const math = require('remark-math');
const katex = require('rehype-katex');


/** @type {import('@docusaurus/types').Config} */
const config = {
  title: 'MontyCarlo',
  tagline: 'The Django of particle simulators.',
  url: 'https://your-docusaurus-test-site.com',
  baseUrl: '/MontyCarlo/',
  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'warn',
  favicon: 'img/favicon.ico',
  organizationName: '', // Usually your GitHub org/user name.
  projectName: 'MontyCarlo', // Usually your repo name.


  stylesheets: [
    {
        href: "https://cdn.jsdelivr.net/npm/katex@0.13.11/dist/katex.min.css",
        integrity: "sha384-Um5gpz1odJg5Z4HAmzPtgZKdTBHZdw8S29IecapCSB31ligYPhHQZMIlWLYQGVoc",
        crossorigin: "anonymous",
    },
],

  plugins: [
        [
         '@docusaurus/plugin-content-docs', 
         
          {
            id: 'physics',
            path: "./physics",
            routeBasePath: 'physics',        
            sidebarPath: require.resolve('./sidebars/physics.js'),
            include: ["**/*.md"],
            remarkPlugins: [math],
            rehypePlugins: [katex],
            // ... other options      
          }, 
        ],  


        [
          '@docusaurus/plugin-content-docs', 
           {
             id: 'tutorial',
             path: "./tutorial",
             routeBasePath: 'tutorial',        
             sidebarPath: require.resolve('./sidebars.js'),
             include: ["**/*.md"],
             // ... other options      
           }, 
         ],  
      ],

  presets: [
    [
      '@docusaurus/preset-classic',
      /** @type {import('@docusaurus/preset-classic').Options} */
      (
        {

          docs: {
            sidebarPath: require.resolve('./sidebars.js'),
            path:'docs',
            // Please change this to your repo.
            editUrl: 'https://github.com/facebook/docusaurus/edit/main/website/',
          },



          blog: {
            showReadingTime: true,
            // Please change this to your repo.
            editUrl:
              'https://github.com/facebook/docusaurus/edit/main/website/blog/',
          },

          theme: {
            customCss: require.resolve('./src/css/custom.css'),
          },
        }
      ),
    ],
  ],



  themeConfig:
    /** @type {import('@docusaurus/preset-classic').ThemeConfig} */
    ({
      navbar: {
        title: 'MontyCarlo.docs',
        logo: {
          alt: 'My Site Logo',
          src: 'img/logo.svg',
        },
        items: [

          {
            to: '/tutorial/intro',
            label: 'Tutorial',
            position: 'left'
          },

          {
            type: 'doc',
            docId: 'intro',
            position: 'left',
            label: 'Documentation',
          },



          {to: '/physics/intro', label: 'Physics', position: 'left'},


          {
            href: 'https://github.com/ruifilipecampos/MontyCarlo',
            label: 'GitHub', 
            position: 'right',
          },

        ],
      },


      footer: {
        style: 'dark',
        copyright: `Copyright © ${new Date().getFullYear()} MontyCarlo.docs, Inc. Built with Docusaurus.`,
      },
      prism: {
        theme: lightCodeTheme,
        darkTheme: darkCodeTheme,
      },
    }),
};

module.exports = config;
