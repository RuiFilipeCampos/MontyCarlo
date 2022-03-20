import React, { useEffect, useState } from 'react'
import Link from 'next/link'
import Layout from '../components/Layout'
import { Text, Flex, VStack, HStack, Spacer, Box } from '@chakra-ui/layout'
import { PlusSquareIcon, ViewIcon, ChevronDownIcon, ChevronRightIcon, MinusIcon, HamburgerIcon } from '@chakra-ui/icons'
import { Button, Icon, IconButton, Input, Textarea, useDisclosure, useForceUpdate } from '@chakra-ui/react'
import Editor from "react-simple-code-editor";

import {
  Menu,
  MenuButton,
  MenuList,
  MenuItem,
  MenuItemOption,
  MenuGroup,
  MenuOptionGroup,
  MenuDivider,
} from '@chakra-ui/react'


import {
  Modal,
  ModalOverlay,
  ModalContent,
  ModalHeader,
  ModalFooter,
  ModalBody,
  ModalCloseButton,
} from '@chakra-ui/react'





const CreateFolderMenuItem = ({ dirChildren, setDirChildren }) => {
  const { isOpen, onOpen, onClose } = useDisclosure()


  let [folderName, setFolderName] = useState('')
  const handleSubmit = () => {
    setDirChildren(
      [...dirChildren, <Directory name={folderName} />]
    )
    onClose()

  }

  return <>
    <MenuItem onClick={onOpen}>New Folder</MenuItem>
    <Modal closeOnOverlayClick={false} isOpen={isOpen} onClose={onClose}>
      <ModalOverlay />
      <ModalContent>
        <ModalHeader>Create a new folder</ModalHeader>
        <ModalCloseButton />
        <ModalBody pb={6}>
          <Input onChange={e => setFolderName(e.target.value)} />
        </ModalBody>

        <ModalFooter>
          <Button onClick={handleSubmit} colorScheme='blue' mr={3}>
            Save
          </Button>
          <Button onClick={onClose}>Cancel</Button>
        </ModalFooter>
      </ModalContent>
    </Modal>
  </>
}

const FolderMenu = ({ dirChildren, setDirChildren }) => {



  return <Menu>
    <MenuButton
      as={IconButton}
      aria-label='Options'
      icon={<HamburgerIcon />}
      variant='outline'
    />
    <MenuList>
      <MenuItem>New File</MenuItem>
      <CreateFolderMenuItem
        dirChildren={dirChildren}
        setDirChildren={setDirChildren}
      />
      <MenuItem>Delete</MenuItem>
    </MenuList>
  </Menu>

}






const NewFolderButton = () => {
  let [folderName, setFolderName] = useState('')

  const onSave = () => {


  }

  return <PopUp buttonText="New Folder">
    <Input onChange={e => setFolderName(e.target.value)} />
  </PopUp>
}








































const get = (endpoint: string) => {
  return fetch(endpoint).then(response => response.json())
}

const delete_request = (endpoint: string) => {
  return fetch(endpoint, {method:'DELETE'})
}



const patch = (endpoint: string, payload) => {
  return fetch(endpoint,
    {
      method: 'PATCH',
      body: JSON.stringify(payload),
      headers: {
        'Content-type': 'application/json; charset=UTF-8',
      }
    }
  ).then(response => response.json())
}

const post = (endpoint: string, payload) => {
  return fetch(endpoint,
    {
      method: 'POST',
      body: JSON.stringify(payload),
      headers: {
        'Content-type': 'application/json; charset=UTF-8',
      }
    }
  ).then(response => response.json())
}



const get_project = (id: number) => get(`/api/projects/${id}/`)

const get_folder = (project_id: number, folder_id: number) => {
  return get(`/api/projects/${project_id}/folders/${folder_id}/`)
}

const delete_folder = (project_id: number, folder_id: number) => {
  return delete_request(`/api/projects/${project_id}/folders/${folder_id}/`)
}


const patch_folder = (project_id: number, folder_id: number, payload) => {
  return patch(`/api/projects/${project_id}/folders/${folder_id}/`, payload)
}




const openFile = (project_id, file_id) => { }

const File = ({ file_name, project_id, file_id }) => {
  const [state, setState] = useState({
    file_name: file_name,
    project_id: project_id,
    file_id: file_id
  })

  return <Flex px={2}><Text textColor="MenuText" css={{ cursor: "pointer" }}
    onClick={(e) => openFile(state.project_id, state.file_id)}
    textAlign="left" variant="link">{state.file_name}
  </Text></Flex>
}



import { ContextMenu } from 'chakra-ui-contextmenu';




function PopUp({ children, buttonText }) {
  const { isOpen, onOpen, onClose } = useDisclosure()

  return (
    <>
      <Button onClick={onOpen}>{buttonText}</Button>
      <Modal closeOnOverlayClick={false} isOpen={isOpen} onClose={onClose}>
        <ModalOverlay />
        <ModalContent>
          <ModalHeader>Create a new folder</ModalHeader>
          <ModalCloseButton />
          <ModalBody pb={6}>
            {children}
          </ModalBody>

          <ModalFooter>
            <Button colorScheme='blue' mr={3}>
              Save
            </Button>
            <Button onClick={onClose}>Cancel</Button>
          </ModalFooter>
        </ModalContent>
      </Modal>
    </>
  )
}




function PopUpMenuItem({ children, buttonText }) {
  const { isOpen, onOpen, onClose } = useDisclosure()

  return (
    <>
      <MenuItem onClick={onOpen}>{buttonText}</MenuItem>
      <Modal closeOnOverlayClick={false} isOpen={isOpen} onClose={onClose}>
        <ModalOverlay />
        <ModalContent>
          <ModalHeader>Create a new folder</ModalHeader>
          <ModalCloseButton />
          <ModalBody pb={6}>
            {children}
          </ModalBody>

          <ModalFooter>
            <Button colorScheme='blue' mr={3}>
              Save
            </Button>
            <Button onClick={onClose}>Cancel</Button>
          </ModalFooter>
        </ModalContent>
      </Modal>
    </>
  )
}




const FolderContextMenu = (
  {children, onClickDeleteFolder, onOpenFile, onOpenFolder, project_id, folder_id}
) => {



  return <ContextMenu renderMenu={() => (
          <MenuList>
            <MenuItem onClick={onOpenFolder}>
              New Folder
            </MenuItem>
            <MenuItem onClick={onOpenFile}>
              New File
            </MenuItem>
            <MenuItem onClick={onClickDeleteFolder}>
              Delete Folder
            </MenuItem>
          </MenuList>
        )}>
          {ref => (
            <div ref={ref}>{children}</div>
          )}
  </ContextMenu>

}




const RenderModalInput = ({ onSubmit, title, isOpen, onClose }) => {
  const [value, setValue] = useState('')

  return <Modal 
      closeOnOverlayClick={false} 
      isOpen={isOpen} 
      onClose={onClose}
    >
  <ModalOverlay />
  <ModalContent>
    <ModalHeader>{title}</ModalHeader>
    <ModalCloseButton />
    <ModalBody pb={6}>
        <Input onChange={e=>setValue(e.target.value)}/>
    </ModalBody>

    <ModalFooter>
      <Button colorScheme='blue' onClick={(e)=>onSubmit(value)} mr={3}>
        Save
      </Button>
      <Button onClick={onClose}>Cancel</Button>
    </ModalFooter>
  </ModalContent>
  </Modal>

  
}


const post_folder = (name, project_id, folder_id)=>{
  return post(
    `/api/projects/${project_id}/folders/`,
    {
      name:name,
      parent:folder_id,
      project:project_id,
      is_open:false
    }
  )

}

const Directory = ({ folder_id, project_id, forceUpdateParent, setforceUpdateParent }) => {


  const [state, setState] = useState({
    name: '',
    is_open: false,
    directories: [],
    code_files: []
  })

  const [forceUpdate, setforceUpdate] = useState(true)
  const [highlight, setHighlight] = useState('black')


  useEffect(() => {
    get_folder(project_id, folder_id).then(payload => setState(payload))
  }, [forceUpdate])



  const change = () => {
    patch_folder(
      project_id,
      folder_id,
      { is_open: !state.is_open }
    ).then(payload => setState(payload)).catch(e=>console.log(e))
  }



  let folderDisclosure = useDisclosure()
  let fileDisclosure = useDisclosure()

  const onClickNewFolder = (name) => {
    post_folder(name, project_id, folder_id).then(
      (e)=>setforceUpdate(!forceUpdate)
    )
    folderDisclosure.onClose()
  }

  const onClickNewFile = (name) => {


    fileDisclosure.onClose()
  }


  const onClickDeleteFolder = () => {
    delete_folder(project_id, folder_id).then(
      (e)=>setforceUpdateParent(!forceUpdateParent())
    ).catch((e)=>console.log(e))

    folderDisclosure.onClose()
  }

  if (state.is_open) return <VStack 
      spacing={.01}
      marginLeft={1} align="left" w="full" 
      onMouseEnter={() => setHighlight('red')} 
      onMouseLeave={() => setHighlight('black')}>

    <RenderModalInput
      title    = "Create new folder"
      isOpen   = {folderDisclosure.isOpen}
      onClose  = {folderDisclosure.onClose} 
      onSubmit = {onClickNewFolder}
      />

    <RenderModalInput 
      title="Create new file"
      isOpen={fileDisclosure.isOpen}
      onClose={fileDisclosure.onClose}
      onSubmit={onClickNewFile}
    />

    <Flex align="center" w="full" >
      <VStack align="left" spacing={0.01}>
        <FolderContextMenu 
          onOpenFolder={folderDisclosure.onOpen} 
          onOpenFile={fileDisclosure.onOpen}
          project_id={project_id}
          folder_id={folder_id}
          onClickDeleteFolder={onClickDeleteFolder}
        >
            <HStack spacing={1} align="center"  css={{ cursor: "pointer" }}>
              <ChevronDownIcon />
              <Text onClick={change}>{state.name}</Text>
              <Spacer />
            </HStack>
        </FolderContextMenu>
      </VStack>
    </Flex>






    <HStack px={1.5} spacing={1} align="right">
      <Box w="3px" h="full" bg={highlight}><></></Box>
      <VStack align="left" w="full" spacing={1}>
        <Spacer />
        {
          state.directories.map(
            (id) => <Directory 
                        forceUpdateParent={() => forceUpdate}
                        setforceUpdateParent={(s) => setforceUpdate(s)}
                        project_id={3} 
                        folder_id={id}
                    />
          )
        }
        {
          state.code_files.map(
            (file_data) => <File
              project_id={3}
              file_id={file_data.id}
              file_name={file_data.name}
            />
          )
        }
      </VStack>
    </HStack>
  </VStack>


  return <VStack marginLeft={3} align="left" w="full">
    <RenderModalInput
      title="Create new folder"
      isOpen={folderDisclosure.isOpen}
      onClose={folderDisclosure.onClose} 
      onSubmit={onClickNewFolder}
      />

    <RenderModalInput 
      title="Create new file"
      isOpen={fileDisclosure.isOpen}
      onClose={fileDisclosure.onClose}
      onSubmit={onClickNewFile}
    />
      <Flex align="center" onClick={change} w="full">
        <VStack align="left" spacing={0.01}>


        <FolderContextMenu 
          onOpenFolder={folderDisclosure.onOpen} 
          onOpenFile={fileDisclosure.onOpen}
          project_id={project_id}
          folder_id={folder_id}
          onClickDeleteFolder={onClickDeleteFolder}
        >
          <HStack align="center" spacing={1}  css={{ cursor: "pointer" }}>
            <ChevronRightIcon />
            <Text onClick={change} >{state.name}</Text>
            <Spacer />
          </HStack>
          </FolderContextMenu>
        </VStack>
      </Flex>
    </VStack>

}











const DirectoryAndFiles = () => {

  const [state, setState] = useState({
    name: 'Project name',
    directories: []
  })

  useEffect(() => {
    get_project(3)
      .then(payload => setState(payload))
      .catch(error => console.log(error))
    }
    , [])




  return <Box > 
    <VStack align="left" h="full" w="full"  spacing={1}>
      <Text>{state.name}</Text>
      <Spacer />
      <Spacer />
      {
        state.directories.map(
          (id) => <Directory
            project_id={3}
            folder_id={id}
          />
        )
      }
    </VStack>
  </Box>
}




































const OpenTerminal = () => <Flex bg="black" h="40%" rounded="lg" w="95%" overflow="hidden" >
  <VStack p={5} align="left">
    <Spacer />
    <Text textColor="white">(your project name) C:\> ls</Text>
    <Text textColor="white">
      {"CHANGELOG.md STYLEGUIDE.md config do.py manage.py reports"}
    </Text>
    <Text textColor="white">(your project name) C:\></Text>
    <Text textColor="white">(your project name) C:\></Text>
    <Text textColor="white">(your project name) C:\></Text>
  </VStack>
</Flex>



const Terminal = ({ terminalMode }) => {
  if (terminalMode) return <OpenTerminal />

  return <Flex bg="black" h="10%" rounded="lg" w="95%" overflow="hidden" p={5} align="left" >
    <Text textColor="white">(your project name) C:\></Text>
  </Flex>
}


function TextCode() {
  return <VStack verticalAlign="center" alignSelf="center" h="full" w="full">
    <Spacer />
    <Textarea h="95%" w="95%" bg="white"
      borderColor="blackAlpha.300"
      focusBorderColor="blackAlpha.400"
    />
    <Spacer />

  </VStack>
}

const Coder = ({ terminalMode }) => {
  if (terminalMode) return <Flex alignItems="center" h="50%" w="full" >
    <TextCode />
  </Flex>

  return <Flex alignItems="center" h="80%" w="full" >
    <TextCode />
  </Flex>
}


const IndexPage = () => {
  useEffect(() => {
    // add a listener to 'message' channel
    global.ipcRenderer.addListener('message', (_event, args) => {
      alert(args)
    })
  }, [])


  let [terminalMode, setTerminalMode] = useState(false)
  const toggleTerminalMode = () => setTerminalMode(!terminalMode)

  return <Layout title="Home | Next.js + TypeScript + Electron Example">
    <Flex  minW="300px" h="inherit" w="15%">
      <VStack h="full" w="full" p={3} align="left" spacing={1}>
        <DirectoryAndFiles />
      </VStack>
    </Flex>
    <Flex minW="1px" w="1px" bg="blackAlpha.300">
      <></>
    </Flex>
    <Spacer />
    <VStack w="85%" bg="Menu">
      <Coder terminalMode={terminalMode} />
      <Terminal terminalMode={terminalMode} />
      <HStack h="10%" w="95%" >
        <Button colorScheme="teal" w="full" > Compile </Button>
        <Button colorScheme="teal" w="full" > Plot Geometry </Button>
        <Button colorScheme="teal" w="full" > Create New Material </Button>
        <Button colorScheme="teal" w="full" onClick={toggleTerminalMode} > Toggle Terminal </Button>
      </HStack>
    </VStack>
  </Layout>
}
export default IndexPage;


